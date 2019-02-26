% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:17:17
% EndTime: 2019-02-26 22:17:17
% DurationCPUTime: 0.24s
% Computational Cost: add. (426->48), mult. (539->67), div. (0->0), fcn. (460->8), ass. (0->46)
t100 = qJ(2) + qJ(3);
t98 = cos(t100);
t135 = qJ(4) * t98;
t139 = pkin(3) + pkin(4);
t97 = sin(t100);
t99 = qJD(2) + qJD(3);
t150 = t99 * t135 + (-t139 * t99 + qJD(4)) * t97;
t142 = qJD(5) - t99;
t101 = sin(qJ(5));
t104 = cos(qJ(5));
t148 = -t101 * t98 + t104 * t97;
t140 = t142 * t148;
t102 = sin(qJ(2));
t136 = pkin(2) * qJD(2);
t125 = t102 * t136;
t149 = -(pkin(9) + r_i_i_C(3) - pkin(8) - pkin(7)) * qJD(1) - t125 + t150;
t116 = t101 * t97 + t104 * t98;
t84 = t142 * t116;
t103 = sin(qJ(1));
t106 = cos(qJ(1));
t114 = qJD(1) * t148;
t77 = -t103 * t114 - t84 * t106;
t113 = qJD(1) * t116;
t78 = t103 * t113 - t140 * t106;
t145 = t77 * r_i_i_C(1) + t78 * r_i_i_C(2);
t79 = t103 * t84 - t106 * t114;
t80 = t140 * t103 + t106 * t113;
t144 = -t79 * r_i_i_C(1) - t80 * r_i_i_C(2);
t143 = -r_i_i_C(1) * t140 + t84 * r_i_i_C(2);
t138 = pkin(2) * t102;
t137 = t98 * t99;
t132 = qJD(4) * t98;
t131 = qJD(1) * t103;
t130 = qJD(1) * t106;
t124 = qJD(1) * t135;
t123 = t97 * t131;
t122 = t97 * t130;
t119 = t103 * t132 + t106 * t124 - t144;
t118 = t106 * t132 + t139 * t123 - t145;
t115 = -qJ(4) * t97 - t139 * t98;
t111 = t115 * t99;
t110 = -t143 + t150;
t105 = cos(qJ(2));
t109 = qJD(1) * (-t105 * pkin(2) - pkin(1) + t115);
t108 = -t105 * t136 + t111;
t1 = [-t80 * r_i_i_C(1) + t79 * r_i_i_C(2) - t149 * t103 + t106 * t109 (-t135 + t138) * t131 + t108 * t106 + t118, -t103 * t124 + t106 * t111 + t118, t106 * t137 - t123, t145, 0; -t78 * r_i_i_C(1) + t77 * r_i_i_C(2) + t103 * t109 + t149 * t106 (-t139 * t97 - t138) * t130 + t108 * t103 + t119, t103 * t111 - t139 * t122 + t119, t103 * t137 + t122, t144, 0; 0, t110 - t125, t110, t99 * t97, t143, 0;];
JaD_transl  = t1;
