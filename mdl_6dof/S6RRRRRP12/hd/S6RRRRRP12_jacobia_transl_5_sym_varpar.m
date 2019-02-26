% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRP12_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP12_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:46:21
% EndTime: 2019-02-26 22:46:21
% DurationCPUTime: 0.45s
% Computational Cost: add. (535->104), mult. (1483->180), div. (0->0), fcn. (1939->14), ass. (0->69)
t51 = sin(qJ(2));
t54 = cos(qJ(2));
t55 = cos(qJ(1));
t72 = cos(pkin(6));
t65 = t55 * t72;
t80 = sin(qJ(1));
t37 = t80 * t51 - t54 * t65;
t38 = t51 * t65 + t80 * t54;
t50 = sin(qJ(3));
t47 = cos(pkin(7));
t81 = cos(qJ(3));
t68 = t47 * t81;
t45 = sin(pkin(7));
t46 = sin(pkin(6));
t77 = t46 * t55;
t70 = t45 * t77;
t15 = t37 * t68 + t38 * t50 + t81 * t70;
t76 = t47 * t50;
t16 = -t37 * t76 + t38 * t81 - t50 * t70;
t30 = -t37 * t45 + t47 * t77;
t49 = sin(qJ(4));
t53 = cos(qJ(4));
t4 = t16 * t53 - t30 * t49;
t48 = sin(qJ(5));
t52 = cos(qJ(5));
t88 = -t15 * t52 + t4 * t48;
t87 = -t15 * t48 - t4 * t52;
t84 = -t16 * t49 - t30 * t53;
t83 = r_i_i_C(3) + pkin(12);
t82 = pkin(10) * t45;
t79 = t45 * t49;
t78 = t45 * t53;
t75 = t50 * t51;
t74 = t50 * t54;
t73 = t51 * t45;
t71 = t46 * t73;
t69 = t46 * t80;
t67 = t81 * t51;
t66 = t81 * t54;
t64 = t72 * t45;
t63 = t45 * t69;
t62 = t72 * t80;
t61 = t52 * r_i_i_C(1) - t48 * r_i_i_C(2) + pkin(4);
t60 = t48 * r_i_i_C(1) + t52 * r_i_i_C(2) + pkin(11);
t59 = t55 * t51 + t54 * t62;
t58 = t59 * t81;
t57 = -t83 * t49 - t61 * t53 - pkin(3);
t56 = t59 * t45 + t47 * t69;
t39 = -t51 * t62 + t55 * t54;
t36 = -t46 * t54 * t45 + t72 * t47;
t35 = (-t47 * t75 + t66) * t46;
t34 = (t47 * t67 + t74) * t46;
t29 = t50 * t64 + (t47 * t74 + t67) * t46;
t28 = -t81 * t64 + (-t47 * t66 + t75) * t46;
t26 = t35 * t53 + t49 * t71;
t24 = -t39 * t76 - t58;
t23 = t39 * t68 - t59 * t50;
t22 = -t37 * t81 - t38 * t76;
t21 = -t37 * t50 + t38 * t68;
t20 = t39 * t81 + (-t59 * t47 + t63) * t50;
t19 = t39 * t50 + t47 * t58 - t81 * t63;
t14 = t29 * t53 + t36 * t49;
t12 = t24 * t53 + t39 * t79;
t10 = t22 * t53 + t38 * t79;
t8 = t20 * t53 + t56 * t49;
t7 = t20 * t49 - t56 * t53;
t2 = t19 * t48 + t8 * t52;
t1 = t19 * t52 - t8 * t48;
t3 = [-t80 * pkin(1) - t38 * pkin(2) - t16 * pkin(3) - t4 * pkin(4) + pkin(9) * t77 + t30 * pkin(10) - t15 * pkin(11) + t87 * r_i_i_C(1) + t88 * r_i_i_C(2) + t83 * t84 (t12 * t52 + t23 * t48) * r_i_i_C(1) + (-t12 * t48 + t23 * t52) * r_i_i_C(2) + t12 * pkin(4) + t24 * pkin(3) + t23 * pkin(11) - t59 * pkin(2) + t39 * t82 + t83 * (t24 * t49 - t39 * t78) t57 * t19 + t60 * t20, -t61 * t7 + t83 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t55 * pkin(1) + t39 * pkin(2) + t20 * pkin(3) + t8 * pkin(4) + pkin(9) * t69 + t56 * pkin(10) + t19 * pkin(11) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t83 * t7 (t10 * t52 + t21 * t48) * r_i_i_C(1) + (-t10 * t48 + t21 * t52) * r_i_i_C(2) + t10 * pkin(4) + t22 * pkin(3) + t21 * pkin(11) - t37 * pkin(2) + t38 * t82 + t83 * (t22 * t49 - t38 * t78) t57 * t15 + t60 * t16, t83 * t4 + t61 * t84, -t88 * r_i_i_C(1) + t87 * r_i_i_C(2), 0; 0 (t26 * t52 + t34 * t48) * r_i_i_C(1) + (-t26 * t48 + t34 * t52) * r_i_i_C(2) + t26 * pkin(4) + t35 * pkin(3) + t34 * pkin(11) + (t54 * pkin(2) + pkin(10) * t73) * t46 + t83 * (t35 * t49 - t53 * t71) t57 * t28 + t60 * t29, t83 * t14 + t61 * (-t29 * t49 + t36 * t53) (-t14 * t48 + t28 * t52) * r_i_i_C(1) + (-t14 * t52 - t28 * t48) * r_i_i_C(2), 0;];
Ja_transl  = t3;
