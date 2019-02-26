% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR6_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR6_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR6_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:56:58
% EndTime: 2019-02-26 21:56:58
% DurationCPUTime: 0.25s
% Computational Cost: add. (361->53), mult. (612->76), div. (0->0), fcn. (520->8), ass. (0->44)
t85 = qJ(4) + qJ(5);
t82 = sin(t85);
t83 = cos(t85);
t87 = sin(qJ(2));
t90 = cos(qJ(2));
t106 = t82 * t90 - t83 * t87;
t86 = sin(qJ(4));
t107 = pkin(4) * t86 + qJ(3);
t89 = cos(qJ(4));
t124 = t89 * pkin(4) + pkin(2) + pkin(3);
t97 = t107 * t90 - t124 * t87;
t128 = t97 * qJD(2) + t87 * qJD(3);
t105 = t82 * t87 + t83 * t90;
t84 = qJD(4) + qJD(5);
t127 = qJD(2) - t84;
t76 = t127 * t105;
t102 = qJD(1) * t106;
t88 = sin(qJ(1));
t91 = cos(qJ(1));
t69 = t88 * t102 + t76 * t91;
t101 = qJD(1) * t105;
t114 = qJD(2) * t90;
t108 = t91 * t114;
t115 = qJD(2) * t87;
t125 = t106 * t84 + t83 * t115;
t70 = t88 * t101 - t82 * t108 + t125 * t91;
t121 = t69 * r_i_i_C(1) + t70 * r_i_i_C(2);
t71 = t91 * t102 - t88 * t76;
t72 = t91 * t101 + (t82 * t114 - t125) * t88;
t120 = -t71 * r_i_i_C(1) - t72 * r_i_i_C(2);
t119 = -t127 * t106 * r_i_i_C(1) - t76 * r_i_i_C(2);
t126 = (pkin(7) - r_i_i_C(3) - pkin(9) - pkin(8)) * qJD(1);
t118 = pkin(4) * qJD(4);
t117 = qJD(1) * t88;
t116 = qJD(1) * t91;
t104 = t86 * t90 - t87 * t89;
t103 = t86 * t87 + t89 * t90;
t100 = t104 * qJD(4);
t99 = -t107 * t87 - t124 * t90;
t96 = qJD(1) * (-pkin(1) + t99);
t95 = (qJD(2) - qJD(4)) * t103;
t94 = t99 * qJD(2) + qJD(3) * t90 + t103 * t118;
t93 = -t104 * t118 + t128;
t1 = [-t72 * r_i_i_C(1) + t71 * r_i_i_C(2) + t91 * t96 + (pkin(4) * t100 - t126 - t128) * t88, -t117 * t97 + t94 * t91 - t121, -t87 * t117 + t108 (t104 * t117 + t95 * t91) * pkin(4) + t121, t121, 0; -t70 * r_i_i_C(1) + t69 * r_i_i_C(2) + t88 * t96 + (t93 + t126) * t91, t97 * t116 + t94 * t88 - t120, t88 * t114 + t87 * t116 (-t104 * t116 + t95 * t88) * pkin(4) + t120, t120, 0; 0, t93 - t119, t115 (-t104 * qJD(2) + t100) * pkin(4) + t119, t119, 0;];
JaD_transl  = t1;
