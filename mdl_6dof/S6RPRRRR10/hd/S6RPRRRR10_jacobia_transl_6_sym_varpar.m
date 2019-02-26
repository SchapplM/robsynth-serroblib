% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRR10_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR10_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobia_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:19:49
% EndTime: 2019-02-26 21:19:50
% DurationCPUTime: 0.35s
% Computational Cost: add. (660->74), mult. (1502->121), div. (0->0), fcn. (1972->16), ass. (0->60)
t52 = cos(qJ(1));
t72 = sin(pkin(13));
t76 = cos(pkin(6));
t64 = t76 * t72;
t74 = cos(pkin(13));
t78 = sin(qJ(1));
t35 = t52 * t64 + t78 * t74;
t49 = sin(qJ(3));
t79 = cos(qJ(3));
t66 = t76 * t74;
t34 = -t52 * t66 + t78 * t72;
t46 = sin(pkin(6));
t73 = sin(pkin(7));
t68 = t46 * t73;
t75 = cos(pkin(7));
t87 = t34 * t75 + t52 * t68;
t21 = t35 * t79 - t87 * t49;
t69 = t46 * t75;
t30 = t34 * t73 - t52 * t69;
t45 = qJ(4) + qJ(5);
t43 = sin(t45);
t44 = cos(t45);
t10 = t21 * t44 + t30 * t43;
t20 = t35 * t49 + t87 * t79;
t47 = sin(qJ(6));
t50 = cos(qJ(6));
t92 = t10 * t47 - t20 * t50;
t91 = -t10 * t50 - t20 * t47;
t48 = sin(qJ(4));
t90 = pkin(4) * t48 + pkin(9);
t83 = r_i_i_C(3) + pkin(12);
t86 = t50 * r_i_i_C(1) - t47 * r_i_i_C(2) + pkin(5);
t85 = -t21 * t43 + t30 * t44;
t56 = t52 * t72 + t78 * t66;
t84 = t56 * t75 - t78 * t68;
t77 = t52 * t46;
t71 = t78 * t46;
t65 = t76 * t73;
t63 = t75 * t74;
t53 = -pkin(11) - pkin(10);
t62 = t47 * r_i_i_C(1) + t50 * r_i_i_C(2) - t53;
t60 = t83 * t10 + t86 * t85;
t36 = t52 * t74 - t78 * t64;
t25 = t36 * t79 - t84 * t49;
t54 = t56 * t73 + t78 * t69;
t13 = t25 * t43 - t54 * t44;
t14 = t25 * t44 + t54 * t43;
t59 = -t86 * t13 + t83 * t14;
t28 = t49 * t65 + (t49 * t63 + t79 * t72) * t46;
t33 = -t74 * t68 + t76 * t75;
t19 = t28 * t44 + t33 * t43;
t58 = t83 * t19 + t86 * (-t28 * t43 + t33 * t44);
t51 = cos(qJ(4));
t42 = t51 * pkin(4) + pkin(3);
t57 = -t83 * t43 - t86 * t44 - t42;
t27 = -t79 * t65 + (t49 * t72 - t63 * t79) * t46;
t24 = t36 * t49 + t84 * t79;
t2 = t14 * t50 + t24 * t47;
t1 = -t14 * t47 + t24 * t50;
t3 = [-t78 * pkin(1) - t35 * pkin(2) - t10 * pkin(5) + t91 * r_i_i_C(1) + t92 * r_i_i_C(2) + qJ(2) * t77 + t20 * t53 - t21 * t42 - t90 * t30 + t83 * t85, t71, t57 * t24 + t62 * t25 (-t25 * t48 + t54 * t51) * pkin(4) + t59, t59, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t52 * pkin(1) + t36 * pkin(2) + t14 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + qJ(2) * t71 + t83 * t13 - t24 * t53 + t25 * t42 + t54 * t90, -t77, t57 * t20 + t62 * t21 (-t21 * t48 + t30 * t51) * pkin(4) + t60, t60, -t92 * r_i_i_C(1) + t91 * r_i_i_C(2); 0, t76, t57 * t27 + t62 * t28 (-t28 * t48 + t33 * t51) * pkin(4) + t58, t58 (-t19 * t47 + t27 * t50) * r_i_i_C(1) + (-t19 * t50 - t27 * t47) * r_i_i_C(2);];
Ja_transl  = t3;
