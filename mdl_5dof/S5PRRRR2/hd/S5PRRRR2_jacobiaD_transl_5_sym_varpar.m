% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a4]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-06-03 15:11
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRRR2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_jacobiaD_transl_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_jacobiaD_transl_5_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5PRRRR2_jacobiaD_transl_5_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-06-03 15:11:34
% EndTime: 2019-06-03 15:11:35
% DurationCPUTime: 0.17s
% Computational Cost: add. (192->49), mult. (310->87), div. (0->0), fcn. (243->8), ass. (0->47)
t95 = cos(qJ(5));
t121 = qJD(5) * t95;
t91 = qJ(3) + qJ(4);
t88 = sin(t91);
t89 = cos(t91);
t90 = qJD(3) + qJD(4);
t92 = sin(qJ(5));
t134 = t89 * t90 * t92 + t88 * t121;
t122 = qJD(5) * t92;
t113 = t88 * t122;
t133 = r_i_i_C(1) * t113 + t134 * r_i_i_C(2);
t93 = sin(qJ(3));
t132 = pkin(1) * t93;
t131 = r_i_i_C(1) * t95;
t130 = r_i_i_C(3) * t89;
t94 = sin(qJ(2));
t129 = t90 * t94;
t128 = t90 * t95;
t97 = cos(qJ(2));
t127 = t90 * t97;
t126 = t95 * t97;
t125 = qJD(2) * t94;
t124 = qJD(2) * t97;
t123 = qJD(3) * t93;
t120 = r_i_i_C(2) * t88 * t92;
t119 = t88 * t129;
t118 = t88 * t128;
t117 = t88 * t127;
t115 = t89 * t127;
t114 = t88 * t125;
t110 = t114 * t131 + t133 * t97;
t108 = qJD(2) * t120;
t106 = qJD(5) * t89 - qJD(2);
t105 = qJD(2) * t89 - qJD(5);
t104 = t97 * t108 + t124 * t130 + t133 * t94;
t103 = t106 * t92;
t102 = -t120 - t130;
t101 = -t88 * t124 - t89 * t129;
t100 = t105 * t94 + t117;
t99 = t89 * r_i_i_C(2) * t121 + t102 * t90 + (t122 * t89 + t118) * r_i_i_C(1);
t96 = cos(qJ(3));
t98 = -pkin(1) * qJD(3) * t96 + (-r_i_i_C(3) * t88 - t89 * t131) * t90;
t75 = -t105 * t126 + (t103 + t118) * t94;
t74 = t106 * t95 * t94 + (t105 * t97 - t119) * t92;
t73 = t100 * t95 + t97 * t103;
t72 = t100 * t92 - t106 * t126;
t1 = [0, t75 * r_i_i_C(1) + t74 * r_i_i_C(2) + t101 * r_i_i_C(3) + (t94 * t123 - t96 * t124) * pkin(1), t98 * t97 + (t102 + t132) * t125 + t110, -t115 * t131 - t94 * t108 + (-t89 * t125 - t117) * r_i_i_C(3) + t110, t72 * r_i_i_C(1) + t73 * r_i_i_C(2); 0, 0, pkin(1) * t123 + t99, t99, (t89 * t128 - t113) * r_i_i_C(2) + t134 * r_i_i_C(1); 0, -t73 * r_i_i_C(1) + t72 * r_i_i_C(2) + (-t114 + t115) * r_i_i_C(3) + (-t97 * t123 - t96 * t125) * pkin(1), (-t88 * t131 - t132) * t124 + t98 * t94 + t104, -r_i_i_C(3) * t119 + t101 * t131 + t104, -t74 * r_i_i_C(1) + t75 * r_i_i_C(2);];
JaD_transl  = t1;
