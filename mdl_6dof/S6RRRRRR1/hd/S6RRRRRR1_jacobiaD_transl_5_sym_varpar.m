% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:47:03
% EndTime: 2019-02-26 22:47:03
% DurationCPUTime: 0.13s
% Computational Cost: add. (398->47), mult. (234->56), div. (0->0), fcn. (150->10), ass. (0->45)
t66 = sin(qJ(2));
t65 = qJ(2) + qJ(3);
t60 = sin(t65);
t64 = qJD(2) + qJD(3);
t59 = qJD(4) + t64;
t55 = qJD(5) + t59;
t62 = qJ(4) + t65;
t58 = qJ(5) + t62;
t54 = cos(t58);
t92 = r_i_i_C(2) * t54;
t53 = sin(t58);
t94 = r_i_i_C(1) * t53;
t75 = t92 + t94;
t73 = t75 * t55;
t56 = sin(t62);
t95 = pkin(4) * t56;
t71 = -t59 * t95 - t73;
t96 = pkin(3) * t64;
t70 = -t60 * t96 + t71;
t88 = pkin(2) * qJD(2);
t98 = -t66 * t88 + t70;
t90 = t54 * t55;
t82 = r_i_i_C(1) * t90;
t57 = cos(t62);
t89 = t57 * t59;
t97 = -pkin(4) * t89 - t82;
t61 = cos(t65);
t76 = -t61 * t96 + t97;
t93 = r_i_i_C(2) * t53;
t91 = r_i_i_C(3) + pkin(10) + pkin(9) + pkin(8) + pkin(7);
t67 = sin(qJ(1));
t87 = qJD(1) * t67;
t69 = cos(qJ(1));
t86 = qJD(1) * t69;
t81 = t55 * t93;
t79 = qJD(1) * t92;
t80 = t67 * t79 + t69 * t81 + t87 * t94;
t68 = cos(qJ(2));
t77 = -t68 * t88 + t76;
t50 = -pkin(3) * t60 - t95;
t74 = -t68 * pkin(2) - pkin(3) * t61 - pkin(4) * t57 - r_i_i_C(1) * t54 - pkin(1) + t93;
t72 = -t69 * t82 + t80;
t48 = t67 * t81;
t47 = -t66 * pkin(2) + t50;
t1 = [-t98 * t67 + (-t91 * t67 + t74 * t69) * qJD(1), -t47 * t87 + t77 * t69 + t80, -t50 * t87 + t76 * t69 + t80 (t56 * t87 - t69 * t89) * pkin(4) + t72, t72, 0; t98 * t69 + (t74 * t67 + t91 * t69) * qJD(1), t48 + t77 * t67 + (t47 - t75) * t86, t48 + t76 * t67 + (t50 - t75) * t86, t48 + t97 * t67 + (-t75 - t95) * t86, -t69 * t79 + t48 + (-t53 * t86 - t67 * t90) * r_i_i_C(1), 0; 0, t98, t70, t71, -t73, 0;];
JaD_transl  = t1;
