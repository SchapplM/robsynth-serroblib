% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR9_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR9_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:35:13
% EndTime: 2019-02-26 22:35:13
% DurationCPUTime: 0.20s
% Computational Cost: add. (290->51), mult. (462->80), div. (0->0), fcn. (578->12), ass. (0->38)
t50 = r_i_i_C(3) + qJ(5);
t30 = sin(pkin(12));
t32 = cos(pkin(12));
t57 = r_i_i_C(1) * t32 - r_i_i_C(2) * t30 + pkin(4);
t36 = cos(qJ(3));
t26 = t36 * pkin(3) + pkin(2);
t29 = qJ(3) + qJ(4);
t27 = sin(t29);
t28 = cos(t29);
t56 = t50 * t27 + t57 * t28 + t26;
t31 = sin(pkin(6));
t34 = sin(qJ(2));
t53 = t31 * t34;
t35 = sin(qJ(1));
t52 = t31 * t35;
t38 = cos(qJ(1));
t51 = t31 * t38;
t49 = cos(pkin(6));
t37 = cos(qJ(2));
t47 = t38 * t49;
t19 = t34 * t47 + t35 * t37;
t8 = t19 * t28 - t27 * t51;
t48 = t35 * t49;
t33 = sin(qJ(3));
t46 = t31 * (pkin(3) * t33 + pkin(8));
t39 = -pkin(10) - pkin(9);
t44 = t30 * r_i_i_C(1) + t32 * r_i_i_C(2) - t39;
t7 = t19 * t27 + t28 * t51;
t43 = t50 * t8 - t57 * t7;
t21 = -t34 * t48 + t38 * t37;
t11 = t21 * t27 - t28 * t52;
t12 = t21 * t28 + t27 * t52;
t42 = -t57 * t11 + t50 * t12;
t16 = t27 * t53 - t49 * t28;
t41 = t50 * (t49 * t27 + t28 * t53) - t57 * t16;
t20 = t38 * t34 + t37 * t48;
t18 = t35 * t34 - t37 * t47;
t1 = [(-t18 * t30 - t32 * t8) * r_i_i_C(1) + (-t18 * t32 + t30 * t8) * r_i_i_C(2) - t8 * pkin(4) - t19 * t26 + t18 * t39 - t35 * pkin(1) - t50 * t7 + t38 * t46, -t20 * t56 + t44 * t21 (-t21 * t33 + t36 * t52) * pkin(3) + t42, t42, t11, 0; (t12 * t32 + t20 * t30) * r_i_i_C(1) + (-t12 * t30 + t20 * t32) * r_i_i_C(2) + t12 * pkin(4) + t21 * t26 - t20 * t39 + t38 * pkin(1) + t35 * t46 + t50 * t11, -t18 * t56 + t44 * t19 (-t19 * t33 - t36 * t51) * pkin(3) + t43, t43, t7, 0; 0 (t44 * t34 + t56 * t37) * t31 (-t33 * t53 + t49 * t36) * pkin(3) + t41, t41, t16, 0;];
Ja_transl  = t1;
