% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP11
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
% Datum: 2019-02-26 22:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRP11_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP11_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:45:38
% EndTime: 2019-02-26 22:45:38
% DurationCPUTime: 0.37s
% Computational Cost: add. (642->92), mult. (1735->161), div. (0->0), fcn. (2265->14), ass. (0->61)
t52 = sin(qJ(2));
t55 = cos(qJ(2));
t56 = cos(qJ(1));
t78 = cos(pkin(6));
t71 = t56 * t78;
t84 = sin(qJ(1));
t37 = t84 * t52 - t55 * t71;
t38 = t52 * t71 + t84 * t55;
t51 = sin(qJ(3));
t77 = cos(pkin(7));
t72 = t51 * t77;
t46 = sin(pkin(7));
t47 = sin(pkin(6));
t80 = t47 * t56;
t75 = t46 * t80;
t85 = cos(qJ(3));
t16 = -t37 * t72 + t38 * t85 - t51 * t75;
t50 = sin(qJ(4));
t54 = cos(qJ(4));
t73 = t47 * t77;
t62 = t37 * t46 - t56 * t73;
t4 = t16 * t54 + t62 * t50;
t3 = t16 * t50 - t62 * t54;
t89 = pkin(5) + r_i_i_C(1);
t67 = t78 * t84;
t61 = t56 * t52 + t55 * t67;
t74 = t47 * t84;
t88 = -t46 * t74 + t61 * t77;
t87 = pkin(10) * t46;
t86 = r_i_i_C(3) + qJ(6) + pkin(12);
t83 = t46 * t50;
t82 = t46 * t54;
t81 = t47 * t55;
t79 = t52 * t46;
t76 = t47 * t79;
t70 = t78 * t46;
t39 = -t52 * t67 + t56 * t55;
t19 = t39 * t51 + t88 * t85;
t49 = sin(qJ(5));
t53 = cos(qJ(5));
t20 = t39 * t85 - t88 * t51;
t57 = t61 * t46 + t84 * t73;
t8 = t20 * t54 + t50 * t57;
t1 = t19 * t53 - t8 * t49;
t66 = t77 * t85;
t45 = t53 * pkin(5) + pkin(4);
t64 = t53 * r_i_i_C(1) - t49 * r_i_i_C(2) + t45;
t63 = t53 * r_i_i_C(2) + t89 * t49 + pkin(11);
t60 = -t46 * t81 + t78 * t77;
t58 = -t86 * t50 - t64 * t54 - pkin(3);
t15 = t37 * t66 + t38 * t51 + t85 * t75;
t35 = (-t52 * t72 + t85 * t55) * t47;
t30 = t51 * t70 + (t85 * t52 + t55 * t72) * t47;
t29 = t47 * t52 * t51 - t66 * t81 - t85 * t70;
t24 = -t39 * t72 - t61 * t85;
t22 = -t37 * t85 - t38 * t72;
t14 = t30 * t54 + t60 * t50;
t13 = t30 * t50 - t60 * t54;
t7 = t20 * t50 - t54 * t57;
t2 = t19 * t49 + t8 * t53;
t5 = [-t84 * pkin(1) - t38 * pkin(2) - t16 * pkin(3) + pkin(9) * t80 - t62 * pkin(10) - t63 * t15 - t86 * t3 - t64 * t4, t24 * pkin(3) - t61 * pkin(2) + t39 * t87 + t64 * (t24 * t54 + t39 * t83) + t63 * (t39 * t66 - t61 * t51) + t86 * (t24 * t50 - t39 * t82) t58 * t19 + t63 * t20, -t64 * t7 + t86 * t8, -t2 * r_i_i_C(2) + t89 * t1, t7; pkin(9) * t74 + t56 * pkin(1) + t39 * pkin(2) + t20 * pkin(3) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t8 * t45 + t86 * t7 + (t49 * pkin(5) + pkin(11)) * t19 + t57 * pkin(10), t38 * t87 - t37 * pkin(2) + t22 * pkin(3) + t86 * (t22 * t50 - t38 * t82) + t64 * (t22 * t54 + t38 * t83) + t63 * (-t37 * t51 + t38 * t66) t58 * t15 + t63 * t16, -t64 * t3 + t86 * t4 (-t15 * t49 - t4 * t53) * r_i_i_C(2) + t89 * (t15 * t53 - t4 * t49) t3; 0, t35 * pkin(3) + t64 * (t35 * t54 + t50 * t76) + t86 * (t35 * t50 - t54 * t76) + (t55 * pkin(2) + pkin(10) * t79 + t63 * (t51 * t55 + t52 * t66)) * t47, t58 * t29 + t63 * t30, -t64 * t13 + t86 * t14 (-t14 * t53 - t29 * t49) * r_i_i_C(2) + t89 * (-t14 * t49 + t29 * t53) t13;];
Ja_transl  = t5;
