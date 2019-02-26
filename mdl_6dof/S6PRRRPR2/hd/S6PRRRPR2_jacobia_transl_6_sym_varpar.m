% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPR2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:11:06
% EndTime: 2019-02-26 20:11:07
% DurationCPUTime: 0.20s
% Computational Cost: add. (310->44), mult. (426->76), div. (0->0), fcn. (532->14), ass. (0->39)
t63 = r_i_i_C(3) + pkin(10) + qJ(5);
t32 = pkin(12) + qJ(6);
t28 = sin(t32);
t29 = cos(t32);
t62 = cos(pkin(12)) * pkin(5) + pkin(4) + t29 * r_i_i_C(1) - t28 * r_i_i_C(2);
t33 = qJ(3) + qJ(4);
t30 = sin(t33);
t31 = cos(t33);
t41 = cos(qJ(3));
t61 = t41 * pkin(3) + t63 * t30 + t62 * t31 + pkin(2);
t35 = sin(pkin(11));
t36 = sin(pkin(6));
t57 = t35 * t36;
t37 = cos(pkin(11));
t56 = t36 * t37;
t40 = sin(qJ(2));
t55 = t36 * t40;
t54 = t36 * t41;
t42 = cos(qJ(2));
t53 = t36 * t42;
t52 = cos(pkin(6));
t51 = t40 * t52;
t50 = t42 * t52;
t20 = t35 * t42 + t37 * t51;
t10 = t20 * t31 - t30 * t56;
t9 = t20 * t30 + t31 * t56;
t48 = t63 * t10 - t62 * t9;
t22 = -t35 * t51 + t37 * t42;
t11 = t22 * t30 - t31 * t57;
t12 = t22 * t31 + t30 * t57;
t47 = -t62 * t11 + t63 * t12;
t17 = t30 * t55 - t52 * t31;
t18 = t52 * t30 + t31 * t55;
t46 = -t62 * t17 + t63 * t18;
t45 = sin(pkin(12)) * pkin(5) + t28 * r_i_i_C(1) + t29 * r_i_i_C(2) + pkin(9) + pkin(8);
t39 = sin(qJ(3));
t21 = t35 * t50 + t37 * t40;
t19 = t35 * t40 - t37 * t50;
t1 = [0, -t21 * t61 + t45 * t22 (-t22 * t39 + t35 * t54) * pkin(3) + t47, t47, t11 (-t12 * t28 + t21 * t29) * r_i_i_C(1) + (-t12 * t29 - t21 * t28) * r_i_i_C(2); 0, -t19 * t61 + t45 * t20 (-t20 * t39 - t37 * t54) * pkin(3) + t48, t48, t9 (-t10 * t28 + t19 * t29) * r_i_i_C(1) + (-t10 * t29 - t19 * t28) * r_i_i_C(2); 1 (t45 * t40 + t61 * t42) * t36 (-t39 * t55 + t52 * t41) * pkin(3) + t46, t46, t17 (-t18 * t28 - t29 * t53) * r_i_i_C(1) + (-t18 * t29 + t28 * t53) * r_i_i_C(2);];
Ja_transl  = t1;
