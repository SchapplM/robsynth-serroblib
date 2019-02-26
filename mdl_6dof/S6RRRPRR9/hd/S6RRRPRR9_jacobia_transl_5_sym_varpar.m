% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR9_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR9_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobia_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:20:27
% EndTime: 2019-02-26 22:20:27
% DurationCPUTime: 0.33s
% Computational Cost: add. (402->90), mult. (1067->155), div. (0->0), fcn. (1393->14), ass. (0->52)
t46 = sin(qJ(2));
t47 = sin(qJ(1));
t50 = cos(qJ(2));
t51 = cos(qJ(1));
t61 = cos(pkin(6));
t57 = t51 * t61;
t30 = t47 * t46 - t50 * t57;
t41 = sin(pkin(7));
t43 = cos(pkin(7));
t42 = sin(pkin(6));
t63 = t42 * t51;
t20 = -t30 * t41 + t43 * t63;
t44 = sin(qJ(5));
t48 = cos(qJ(5));
t40 = sin(pkin(13));
t45 = sin(qJ(3));
t49 = cos(qJ(3));
t60 = cos(pkin(13));
t35 = -t49 * t40 - t45 * t60;
t24 = t35 * t41;
t26 = t35 * t43;
t31 = t46 * t57 + t47 * t50;
t53 = -t45 * t40 + t49 * t60;
t7 = -t24 * t63 - t30 * t26 - t31 * t53;
t70 = t20 * t48 - t7 * t44;
t69 = t20 * t44 + t7 * t48;
t13 = (-t26 * t50 + t46 * t53) * t42 - t61 * t24;
t68 = -r_i_i_C(3) - pkin(11);
t67 = pkin(3) * t45;
t66 = t41 * t44;
t65 = t41 * t48;
t64 = t42 * t47;
t62 = pkin(10) + qJ(4);
t59 = t42 * (t41 * t67 + t62 * t43 + pkin(9));
t58 = t47 * t61;
t54 = t48 * r_i_i_C(1) - t44 * r_i_i_C(2) + pkin(4);
t23 = t53 * t41;
t25 = t53 * t43;
t52 = t23 * t63 + t30 * t25 - t31 * t35;
t32 = -t51 * t46 - t50 * t58;
t33 = -t46 * t58 + t51 * t50;
t10 = -t24 * t64 - t32 * t26 + t33 * t53;
t39 = t49 * pkin(3) + pkin(2);
t29 = -t42 * t50 * t41 + t61 * t43;
t28 = -t62 * t41 + t43 * t67;
t22 = -t32 * t41 + t43 * t64;
t17 = t33 * t26 + t32 * t53;
t15 = t31 * t26 - t30 * t53;
t9 = t23 * t64 + t32 * t25 + t33 * t35;
t2 = t10 * t48 + t22 * t44;
t1 = -t10 * t44 + t22 * t48;
t3 = [-t47 * pkin(1) + t7 * pkin(4) + t69 * r_i_i_C(1) + t70 * r_i_i_C(2) + t30 * t28 - t31 * t39 + t51 * t59 + t68 * t52 (t17 * t48 + t33 * t66) * r_i_i_C(1) + (-t17 * t44 + t33 * t65) * r_i_i_C(2) + t17 * pkin(4) + t32 * t39 - t33 * t28 + t68 * (-t33 * t25 + t32 * t35) -t68 * t10 + t54 * t9 + (-t33 * t45 + (t32 * t43 + t41 * t64) * t49) * pkin(3), t22, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t51 * pkin(1) + t10 * pkin(4) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t32 * t28 + t33 * t39 + t47 * t59 + t68 * t9 (t15 * t48 + t31 * t66) * r_i_i_C(1) + (-t15 * t44 + t31 * t65) * r_i_i_C(2) + t15 * pkin(4) - t30 * t39 - t31 * t28 + t68 * (-t31 * t25 - t30 * t35) t68 * t7 - t54 * t52 + (-t31 * t45 + (-t30 * t43 - t41 * t63) * t49) * pkin(3), -t20, -t70 * r_i_i_C(1) + t69 * r_i_i_C(2), 0; 0 (t54 * (t26 * t46 + t50 * t53) + t50 * t39 + (-t28 + (t44 * r_i_i_C(1) + t48 * r_i_i_C(2)) * t41) * t46 + t68 * (-t25 * t46 + t35 * t50)) * t42, -t68 * t13 + t54 * (t61 * t23 + (t25 * t50 + t35 * t46) * t42) + (t61 * t41 * t49 + (t43 * t49 * t50 - t45 * t46) * t42) * pkin(3), t29 (-t13 * t44 + t29 * t48) * r_i_i_C(1) + (-t13 * t48 - t29 * t44) * r_i_i_C(2), 0;];
Ja_transl  = t3;
