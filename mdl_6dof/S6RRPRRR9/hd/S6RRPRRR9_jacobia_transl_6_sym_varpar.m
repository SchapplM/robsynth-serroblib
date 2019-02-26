% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR9_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR9_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:58:45
% EndTime: 2019-02-26 21:58:46
% DurationCPUTime: 0.21s
% Computational Cost: add. (449->60), mult. (523->93), div. (0->0), fcn. (656->14), ass. (0->44)
t61 = pkin(11) + r_i_i_C(3);
t37 = sin(qJ(6));
t40 = cos(qJ(6));
t65 = t40 * r_i_i_C(1) - t37 * r_i_i_C(2) + pkin(5);
t38 = sin(qJ(2));
t39 = sin(qJ(1));
t41 = cos(qJ(2));
t42 = cos(qJ(1));
t52 = cos(pkin(6));
t49 = t42 * t52;
t21 = t38 * t49 + t39 * t41;
t35 = pkin(12) + qJ(4);
t33 = qJ(5) + t35;
t29 = sin(t33);
t30 = cos(t33);
t36 = sin(pkin(6));
t53 = t36 * t42;
t10 = t21 * t30 - t29 * t53;
t20 = t39 * t38 - t41 * t49;
t64 = t10 * t37 - t20 * t40;
t63 = -t10 * t40 - t20 * t37;
t32 = cos(t35);
t24 = pkin(4) * t32 + cos(pkin(12)) * pkin(3) + pkin(2);
t62 = t61 * t29 + t65 * t30 + t24;
t56 = t36 * t38;
t55 = t36 * t39;
t54 = t36 * t41;
t31 = sin(t35);
t51 = t36 * (pkin(8) + pkin(4) * t31 + sin(pkin(12)) * pkin(3));
t50 = t39 * t52;
t34 = -pkin(10) - pkin(9) - qJ(3);
t47 = t37 * r_i_i_C(1) + t40 * r_i_i_C(2) - t34;
t9 = -t21 * t29 - t30 * t53;
t46 = t61 * t10 + t65 * t9;
t23 = -t38 * t50 + t42 * t41;
t13 = t23 * t29 - t30 * t55;
t14 = t23 * t30 + t29 * t55;
t45 = -t65 * t13 + t61 * t14;
t19 = t52 * t29 + t30 * t56;
t44 = t61 * t19 + t65 * (-t29 * t56 + t52 * t30);
t22 = t42 * t38 + t41 * t50;
t2 = t14 * t40 + t22 * t37;
t1 = -t14 * t37 + t22 * t40;
t3 = [-t39 * pkin(1) - t10 * pkin(5) + t63 * r_i_i_C(1) + t64 * r_i_i_C(2) + t20 * t34 - t21 * t24 + t42 * t51 + t61 * t9, -t22 * t62 + t47 * t23, t22 (-t23 * t31 + t32 * t55) * pkin(4) + t45, t45, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t42 * pkin(1) + t14 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t61 * t13 - t22 * t34 + t23 * t24 + t39 * t51, -t20 * t62 + t47 * t21, t20 (-t21 * t31 - t32 * t53) * pkin(4) + t46, t46, -t64 * r_i_i_C(1) + t63 * r_i_i_C(2); 0 (t47 * t38 + t62 * t41) * t36, -t54 (-t31 * t56 + t52 * t32) * pkin(4) + t44, t44 (-t19 * t37 - t40 * t54) * r_i_i_C(1) + (-t19 * t40 + t37 * t54) * r_i_i_C(2);];
Ja_transl  = t3;
