% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRR4
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRRR4_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR4_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_jacobia_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:20:26
% EndTime: 2019-02-26 20:20:26
% DurationCPUTime: 0.20s
% Computational Cost: add. (287->58), mult. (661->106), div. (0->0), fcn. (849->14), ass. (0->45)
t31 = sin(pkin(13));
t34 = cos(pkin(13));
t39 = sin(qJ(2));
t36 = cos(pkin(6));
t42 = cos(qJ(2));
t52 = t36 * t42;
t22 = -t31 * t39 + t34 * t52;
t32 = sin(pkin(7));
t33 = sin(pkin(6));
t35 = cos(pkin(7));
t56 = t33 * t35;
t17 = -t22 * t32 - t34 * t56;
t30 = qJ(4) + qJ(5);
t28 = sin(t30);
t29 = cos(t30);
t53 = t36 * t39;
t23 = t31 * t42 + t34 * t53;
t38 = sin(qJ(3));
t41 = cos(qJ(3));
t58 = t32 * t33;
t46 = t22 * t35 - t34 * t58;
t8 = t23 * t41 + t46 * t38;
t62 = (t17 * t29 - t8 * t28) * r_i_i_C(1) + (-t17 * t28 - t8 * t29) * r_i_i_C(2);
t25 = -t31 * t53 + t34 * t42;
t24 = -t31 * t52 - t34 * t39;
t45 = t24 * t35 + t31 * t58;
t10 = t25 * t41 + t45 * t38;
t18 = -t24 * t32 + t31 * t56;
t61 = (-t10 * t28 + t18 * t29) * r_i_i_C(1) + (-t10 * t29 - t18 * t28) * r_i_i_C(2);
t49 = t39 * t41;
t50 = t38 * t42;
t57 = t32 * t36;
t16 = t38 * t57 + (t35 * t50 + t49) * t33;
t21 = t36 * t35 - t42 * t58;
t60 = (-t16 * t28 + t21 * t29) * r_i_i_C(1) + (-t16 * t29 - t21 * t28) * r_i_i_C(2);
t59 = r_i_i_C(3) + pkin(11) + pkin(10);
t55 = t35 * t38;
t54 = t35 * t41;
t51 = t38 * t39;
t48 = t41 * t42;
t40 = cos(qJ(4));
t47 = t40 * pkin(4) + t29 * r_i_i_C(1) - t28 * r_i_i_C(2) + pkin(3);
t37 = sin(qJ(4));
t44 = (pkin(4) * t37 + r_i_i_C(1) * t28 + r_i_i_C(2) * t29 + pkin(9)) * t32;
t1 = [0, t24 * pkin(2) + t47 * (t24 * t41 - t25 * t55) + t25 * t44 + t59 * (t24 * t38 + t25 * t54) t59 * t10 + t47 * (-t25 * t38 + t45 * t41) (-t10 * t37 + t18 * t40) * pkin(4) + t61, t61, 0; 0, t22 * pkin(2) + t47 * (t22 * t41 - t23 * t55) + t23 * t44 + t59 * (t22 * t38 + t23 * t54) t59 * t8 + t47 * (-t23 * t38 + t46 * t41) (t17 * t40 - t37 * t8) * pkin(4) + t62, t62, 0; 1 (t59 * (t35 * t49 + t50) + t47 * (-t35 * t51 + t48) + t42 * pkin(2) + t39 * t44) * t33, t59 * t16 + t47 * (t41 * t57 + (t35 * t48 - t51) * t33) (-t16 * t37 + t21 * t40) * pkin(4) + t60, t60, 0;];
Ja_transl  = t1;
