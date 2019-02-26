% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRR6
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRRR6_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR6_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobia_transl_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:22:01
% EndTime: 2019-02-26 20:22:01
% DurationCPUTime: 0.34s
% Computational Cost: add. (640->118), mult. (1862->221), div. (0->0), fcn. (2442->16), ass. (0->83)
t93 = r_i_i_C(3) + pkin(12);
t46 = sin(pkin(8));
t92 = t46 * pkin(11);
t47 = sin(pkin(7));
t91 = t47 * pkin(10);
t48 = sin(pkin(6));
t51 = cos(pkin(7));
t55 = sin(qJ(2));
t58 = cos(qJ(3));
t77 = t55 * t58;
t54 = sin(qJ(3));
t59 = cos(qJ(2));
t78 = t54 * t59;
t37 = (-t51 * t77 - t78) * t48;
t90 = t37 * t46;
t89 = t46 * t47;
t52 = sin(qJ(5));
t88 = t46 * t52;
t56 = cos(qJ(5));
t87 = t46 * t56;
t50 = cos(pkin(8));
t86 = t47 * t50;
t85 = t47 * t55;
t49 = cos(pkin(14));
t84 = t48 * t49;
t53 = sin(qJ(4));
t83 = t50 * t53;
t57 = cos(qJ(4));
t82 = t50 * t57;
t81 = t51 * t54;
t80 = t51 * t58;
t79 = t54 * t55;
t76 = t58 * t59;
t75 = cos(pkin(6));
t74 = sin(pkin(14));
t73 = t47 * t84;
t72 = t48 * t85;
t71 = t48 * t74;
t70 = t49 * t75;
t69 = t75 * t47;
t40 = -t74 * t55 + t59 * t70;
t41 = t55 * t70 + t74 * t59;
t24 = -t41 * t54 + (t40 * t51 - t73) * t58;
t35 = -t40 * t47 - t51 * t84;
t68 = t24 * t50 + t35 * t46;
t65 = t75 * t74;
t43 = t49 * t59 - t55 * t65;
t42 = -t49 * t55 - t59 * t65;
t60 = t42 * t51 + t47 * t71;
t26 = -t43 * t54 + t60 * t58;
t36 = -t42 * t47 + t51 * t71;
t67 = t26 * t50 + t36 * t46;
t33 = t58 * t69 + (t51 * t76 - t79) * t48;
t39 = -t48 * t59 * t47 + t75 * t51;
t66 = t33 * t50 + t39 * t46;
t64 = t56 * r_i_i_C(1) - t52 * r_i_i_C(2) + pkin(4);
t28 = -t40 * t54 - t41 * t80;
t19 = -t28 * t46 + t41 * t86;
t63 = t28 * t50 + t41 * t89;
t30 = -t42 * t54 - t43 * t80;
t20 = -t30 * t46 + t43 * t86;
t62 = t30 * t50 + t43 * t89;
t61 = t37 * t50 + t46 * t72;
t38 = (-t51 * t79 + t76) * t48;
t34 = t54 * t69 + (t51 * t78 + t77) * t48;
t32 = t50 * t72 - t90;
t31 = t42 * t58 - t43 * t81;
t29 = t40 * t58 - t41 * t81;
t27 = t43 * t58 + t60 * t54;
t25 = t40 * t81 + t41 * t58 - t54 * t73;
t23 = -t33 * t46 + t39 * t50;
t22 = t38 * t57 + t61 * t53;
t18 = t33 * t57 - t34 * t83;
t16 = -t26 * t46 + t36 * t50;
t15 = -t24 * t46 + t35 * t50;
t14 = t34 * t57 + t66 * t53;
t12 = t26 * t57 - t27 * t83;
t10 = t24 * t57 - t25 * t83;
t8 = t31 * t57 + t62 * t53;
t6 = t29 * t57 + t63 * t53;
t4 = t27 * t57 + t67 * t53;
t2 = t25 * t57 + t68 * t53;
t1 = [0 (t20 * t52 + t8 * t56) * r_i_i_C(1) + (t20 * t56 - t8 * t52) * r_i_i_C(2) + t8 * pkin(4) + t31 * pkin(3) + t42 * pkin(2) + t43 * t91 + t93 * (t31 * t53 - t62 * t57) + t20 * pkin(11) (t12 * t56 + t27 * t88) * r_i_i_C(1) + (-t12 * t52 + t27 * t87) * r_i_i_C(2) + t12 * pkin(4) + t26 * pkin(3) + t27 * t92 + t93 * (t26 * t53 + t27 * t82) t93 * t4 + t64 * (-t27 * t53 + t67 * t57) (t16 * t56 - t4 * t52) * r_i_i_C(1) + (-t16 * t52 - t4 * t56) * r_i_i_C(2), 0; 0 (t19 * t52 + t6 * t56) * r_i_i_C(1) + (t19 * t56 - t6 * t52) * r_i_i_C(2) + t6 * pkin(4) + t29 * pkin(3) + t40 * pkin(2) + t41 * t91 + t93 * (t29 * t53 - t63 * t57) + t19 * pkin(11) (t10 * t56 + t25 * t88) * r_i_i_C(1) + (-t10 * t52 + t25 * t87) * r_i_i_C(2) + t10 * pkin(4) + t24 * pkin(3) + t25 * t92 + t93 * (t24 * t53 + t25 * t82) t93 * t2 + t64 * (-t25 * t53 + t68 * t57) (t15 * t56 - t2 * t52) * r_i_i_C(1) + (-t15 * t52 - t2 * t56) * r_i_i_C(2), 0; 1 (t22 * t56 + t32 * t52) * r_i_i_C(1) + (-t22 * t52 + t32 * t56) * r_i_i_C(2) + t22 * pkin(4) + t38 * pkin(3) - pkin(11) * t90 + t93 * (t38 * t53 - t61 * t57) + (t59 * pkin(2) + (pkin(11) * t50 + pkin(10)) * t85) * t48 (t18 * t56 + t34 * t88) * r_i_i_C(1) + (-t18 * t52 + t34 * t87) * r_i_i_C(2) + t18 * pkin(4) + t33 * pkin(3) + t34 * t92 + t93 * (t33 * t53 + t34 * t82) t93 * t14 + t64 * (-t34 * t53 + t66 * t57) (-t14 * t52 + t23 * t56) * r_i_i_C(1) + (-t14 * t56 - t23 * t52) * r_i_i_C(2), 0;];
Ja_transl  = t1;
