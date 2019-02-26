% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR7_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR7_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:41:12
% EndTime: 2019-02-26 21:41:13
% DurationCPUTime: 0.22s
% Computational Cost: add. (260->54), mult. (522->76), div. (0->0), fcn. (436->8), ass. (0->43)
t69 = qJ(4) + pkin(10);
t67 = sin(t69);
t68 = cos(t69);
t72 = sin(qJ(2));
t75 = cos(qJ(2));
t90 = t67 * t75 - t68 * t72;
t74 = cos(qJ(4));
t108 = pkin(4) * t74 + pkin(2) + pkin(3);
t71 = sin(qJ(4));
t94 = pkin(4) * t71 + qJ(3);
t82 = t108 * t72 - t94 * t75;
t112 = -t82 * qJD(2) + t72 * qJD(3);
t110 = qJD(2) - qJD(4);
t89 = t67 * t72 + t68 * t75;
t61 = t110 * t89;
t111 = (pkin(7) - r_i_i_C(3) - qJ(5) - pkin(8)) * qJD(1);
t102 = qJD(2) * t72;
t109 = t90 * qJD(4) + t68 * t102;
t105 = pkin(4) * qJD(4);
t73 = sin(qJ(1));
t104 = qJD(1) * t73;
t76 = cos(qJ(1));
t103 = qJD(1) * t76;
t101 = qJD(2) * t75;
t97 = t76 * t101;
t86 = qJD(1) * t90;
t56 = t61 * t76 + t73 * t86;
t85 = qJD(1) * t89;
t57 = t109 * t76 - t67 * t97 + t73 * t85;
t93 = t56 * r_i_i_C(1) + t57 * r_i_i_C(2);
t58 = -t73 * t61 + t76 * t86;
t59 = t76 * t85 + (t67 * t101 - t109) * t73;
t92 = -t58 * r_i_i_C(1) - t59 * r_i_i_C(2);
t91 = -t110 * t90 * r_i_i_C(1) - t61 * r_i_i_C(2);
t88 = t71 * t75 - t72 * t74;
t87 = t71 * t72 + t74 * t75;
t84 = t88 * qJD(4);
t83 = -t108 * t75 - t94 * t72;
t80 = -qJD(5) + (-pkin(1) + t83) * qJD(1);
t79 = t110 * t87;
t78 = t83 * qJD(2) + qJD(3) * t75 + t87 * t105;
t77 = -t88 * t105 + t112;
t1 = [-t59 * r_i_i_C(1) + t58 * r_i_i_C(2) + t80 * t76 + (pkin(4) * t84 - t111 - t112) * t73, t82 * t104 + t78 * t76 - t93, -t72 * t104 + t97 (t88 * t104 + t79 * t76) * pkin(4) + t93, -t103, 0; -t57 * r_i_i_C(1) + t56 * r_i_i_C(2) + t80 * t73 + (t77 + t111) * t76, -t103 * t82 + t78 * t73 - t92, t73 * t101 + t72 * t103 (-t88 * t103 + t79 * t73) * pkin(4) + t92, -t104, 0; 0, t77 - t91, t102 (-t88 * qJD(2) + t84) * pkin(4) + t91, 0, 0;];
JaD_transl  = t1;
