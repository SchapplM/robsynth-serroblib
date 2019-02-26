% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP8
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRP8_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP8_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:43:48
% EndTime: 2019-02-26 22:43:48
% DurationCPUTime: 0.30s
% Computational Cost: add. (478->68), mult. (783->110), div. (0->0), fcn. (994->12), ass. (0->48)
t74 = r_i_i_C(3) + qJ(6);
t55 = sin(qJ(5));
t59 = cos(qJ(5));
t89 = pkin(5) + r_i_i_C(1);
t94 = t55 * t74 + t59 * t89 + pkin(4);
t88 = pkin(11) + r_i_i_C(2);
t57 = sin(qJ(2));
t58 = sin(qJ(1));
t61 = cos(qJ(2));
t62 = cos(qJ(1));
t72 = cos(pkin(6));
t69 = t62 * t72;
t44 = t57 * t69 + t58 * t61;
t53 = qJ(3) + qJ(4);
t51 = sin(t53);
t52 = cos(t53);
t54 = sin(pkin(6));
t77 = t54 * t62;
t26 = t44 * t52 - t51 * t77;
t43 = t57 * t58 - t61 * t69;
t1 = t26 * t55 - t43 * t59;
t91 = t26 * t59 + t43 * t55;
t60 = cos(qJ(3));
t50 = pkin(3) * t60 + pkin(2);
t90 = pkin(4) * t52 + t51 * t88 + t50;
t81 = t52 * t55;
t80 = t52 * t59;
t79 = t54 * t57;
t78 = t54 * t58;
t76 = t55 * t61;
t75 = t59 * t61;
t70 = t58 * t72;
t56 = sin(qJ(3));
t68 = t54 * (pkin(3) * t56 + pkin(8));
t25 = -t44 * t51 - t52 * t77;
t66 = t25 * t94 + t88 * t26;
t46 = -t57 * t70 + t61 * t62;
t29 = t46 * t51 - t52 * t78;
t30 = t46 * t52 + t51 * t78;
t65 = -t29 * t94 + t88 * t30;
t40 = t51 * t72 + t52 * t79;
t64 = t88 * t40 + t94 * (-t51 * t79 + t52 * t72);
t63 = -pkin(10) - pkin(9);
t45 = t57 * t62 + t61 * t70;
t23 = t40 * t55 + t54 * t75;
t6 = t30 * t59 + t45 * t55;
t5 = t30 * t55 - t45 * t59;
t2 = [-t58 * pkin(1) - pkin(4) * t26 - t1 * t74 + t88 * t25 + t43 * t63 - t44 * t50 + t62 * t68 - t89 * t91, -t46 * t63 + t74 * (-t45 * t81 - t46 * t59) + t89 * (-t45 * t80 + t46 * t55) - t90 * t45 (-t46 * t56 + t60 * t78) * pkin(3) + t65, t65, -t5 * t89 + t6 * t74, t5; t62 * pkin(1) + t30 * pkin(4) + t29 * t88 - t45 * t63 + t46 * t50 + t5 * t74 + t58 * t68 + t6 * t89, -t44 * t63 + t89 * (-t43 * t80 + t44 * t55) + t74 * (-t43 * t81 - t44 * t59) - t90 * t43 (-t44 * t56 - t60 * t77) * pkin(3) + t66, t66, -t1 * t89 + t74 * t91, t1; 0 (t89 * (t52 * t75 + t55 * t57) + t74 * (t52 * t76 - t57 * t59) - t57 * t63 + t90 * t61) * t54 (-t56 * t79 + t60 * t72) * pkin(3) + t64, t64, t74 * (t40 * t59 - t54 * t76) - t89 * t23, t23;];
Ja_transl  = t2;
