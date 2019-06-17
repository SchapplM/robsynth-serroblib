% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5RPRRR1
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
%   pkin=[dummy]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-06-12 14:37
% Revision: aab8d7cd0cba739f5e0ec8d53b8419901d1154b0 (2019-06-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRRR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_jacobiaD_transl_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_jacobiaD_transl_5_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_jacobiaD_transl_5_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-06-12 14:37:34
% EndTime: 2019-06-12 14:37:34
% DurationCPUTime: 0.28s
% Computational Cost: add. (185->75), mult. (600->142), div. (0->0), fcn. (566->8), ass. (0->56)
t102 = cos(qJ(4));
t146 = -qJD(5) * t102 + qJD(3);
t101 = cos(qJ(5));
t98 = sin(qJ(4));
t137 = qJD(4) * t98;
t97 = sin(qJ(5));
t109 = t146 * t101 + t97 * t137;
t110 = t101 * t137 - t146 * t97;
t127 = qJD(4) * t102;
t145 = t110 * r_i_i_C(1) - t109 * r_i_i_C(2) - r_i_i_C(3) * t127;
t103 = cos(qJ(3));
t99 = sin(qJ(3));
t140 = t102 * t99;
t144 = (t101 * t140 - t103 * t97) * r_i_i_C(1) - (t101 * t103 + t97 * t140) * r_i_i_C(2) + t99 * t98 * r_i_i_C(3);
t143 = t97 * r_i_i_C(1);
t142 = t97 * r_i_i_C(2);
t100 = sin(qJ(1));
t141 = t100 * t98;
t104 = cos(qJ(1));
t139 = t104 * t98;
t138 = qJD(3) * t98;
t136 = qJD(5) * t97;
t135 = qJD(5) * t99;
t134 = t100 * t103;
t133 = t104 * t102;
t132 = qJD(1) * t100;
t131 = qJD(1) * t104;
t130 = qJD(3) * t102;
t129 = qJD(3) * t103;
t128 = qJD(3) * t104;
t126 = qJD(5) * t101;
t122 = t98 * t134;
t121 = t99 * t131;
t120 = -qJD(4) * t103 + qJD(1);
t119 = qJD(1) * t103 - qJD(4);
t118 = -qJD(5) + t130;
t117 = -t101 * r_i_i_C(1) + t142;
t92 = t103 * t133 + t141;
t88 = -t100 * t99 * t130 + t92 * qJD(1) - qJD(4) * t122 - t104 * t127;
t116 = -t100 * t135 - t88;
t86 = t120 * t139 + (-t119 * t100 - t99 * t128) * t102;
t115 = t104 * t135 + t86;
t114 = t101 * t118;
t111 = t120 * t102 + t99 * t138;
t90 = t102 * t134 - t139;
t108 = -qJD(5) * t90 + t100 * t129 + t121;
t107 = -qJD(5) * t92 + t103 * t128 - t99 * t132;
t106 = -r_i_i_C(1) * t114 - r_i_i_C(3) * t138 + t118 * t142;
t105 = t106 * t103 + t145 * t99;
t91 = t100 * t102 - t103 * t139;
t89 = -t122 - t133;
t87 = t111 * t100 - t119 * t139;
t85 = t111 * t104 + t119 * t141;
t84 = t115 * t101 + t107 * t97;
t83 = t107 * t101 - t115 * t97;
t1 = [(-t88 * t101 - t97 * t121 + t90 * t136) * r_i_i_C(1) + (-t101 * t121 + t90 * t126 + t88 * t97) * r_i_i_C(2) + t87 * r_i_i_C(3) + t104 * qJD(2) + ((-t99 * t126 - t97 * t129) * r_i_i_C(1) + (-t101 * t129 + t97 * t135) * r_i_i_C(2) - qJD(1) * qJ(2)) * t100, t131, t105 * t104 + t144 * t132, t86 * r_i_i_C(3) + (-t91 * t126 - t85 * t97) * r_i_i_C(2) + (t101 * t85 - t91 * t136) * r_i_i_C(1), t83 * r_i_i_C(1) - t84 * r_i_i_C(2); t84 * r_i_i_C(1) + t83 * r_i_i_C(2) - t85 * r_i_i_C(3) + qJ(2) * t131 + t100 * qJD(2), t132, t105 * t100 - t144 * t131, t88 * r_i_i_C(3) + (-t89 * t126 - t87 * t97) * r_i_i_C(2) + (t101 * t87 - t89 * t136) * r_i_i_C(1), (t116 * r_i_i_C(1) - t108 * r_i_i_C(2)) * t97 + (t108 * r_i_i_C(1) + t116 * r_i_i_C(2)) * t101; 0, 0, -t145 * t103 + t106 * t99, (t117 * t99 * qJD(4) + r_i_i_C(3) * t129) * t102 + (t117 * t129 + (-qJD(4) * r_i_i_C(3) + (t101 * r_i_i_C(2) + t143) * qJD(5)) * t99) * t98, (-r_i_i_C(2) * t114 - t118 * t143) * t103 + (t109 * r_i_i_C(1) + t110 * r_i_i_C(2)) * t99;];
JaD_transl  = t1;
