% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S4RRRP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 13:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4RRRP3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRRP3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_jacobiaD_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:54:37
	% EndTime: 2019-12-29 13:54:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:54:37
	% EndTime: 2019-12-29 13:54:37
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:54:30
	% EndTime: 2019-12-29 13:54:30
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (74->0), mult. (74->0), div. (30->0), fcn. (44->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:54:37
	% EndTime: 2019-12-29 13:54:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:54:37
	% EndTime: 2019-12-29 13:54:38
	% DurationCPUTime: 1.38s
	% Computational Cost: add. (2232->71), mult. (2645->160), div. (674->13), fcn. (3179->7), ass. (0->78)
	t109 = qJ(1) + qJ(2);
	t101 = sin(t109);
	t102 = cos(t109);
	t98 = 0.1e1 / t102;
	t150 = t101 * t98;
	t138 = qJD(3) * t101;
	t111 = cos(qJ(3));
	t106 = 0.1e1 / t111;
	t110 = sin(qJ(3));
	t105 = t110 ^ 2;
	t107 = 0.1e1 / t111 ^ 2;
	t141 = t105 * t107;
	t127 = t101 * t141;
	t137 = qJD(3) * t111;
	t103 = qJD(1) + qJD(2);
	t142 = t103 * t110;
	t96 = t101 ^ 2;
	t94 = t96 * t141 + 0.1e1;
	t92 = 0.1e1 / t94;
	t75 = (-(-t101 * t137 - t102 * t142) * t106 + qJD(3) * t127) * t92;
	t162 = t75 - t138;
	t97 = t102 ^ 2;
	t161 = t103 * (0.1e1 / t97 * t96 + 0.1e1) * t150;
	t144 = t101 * t110;
	t91 = atan2(-t144, -t111);
	t86 = sin(t91);
	t131 = t86 * t144;
	t87 = cos(t91);
	t84 = -t111 * t87 - t131;
	t81 = 0.1e1 / t84;
	t82 = 0.1e1 / t84 ^ 2;
	t159 = 0.2e1 * t110;
	t158 = t92 - 0.1e1;
	t128 = t82 * t137;
	t145 = t101 * t103;
	t133 = t82 * t145;
	t143 = t102 * t103;
	t148 = t111 * t86;
	t151 = t101 * t75;
	t70 = t162 * t148 + (-t86 * t143 + (qJD(3) - t151) * t87) * t110;
	t156 = t70 * t81 * t82;
	t78 = t105 * t82 * t97 + 0.1e1;
	t157 = (t97 * t110 * t128 + (-t102 * t133 - t97 * t156) * t105) / t78 ^ 2;
	t76 = 0.1e1 / t78;
	t155 = t76 * t82;
	t104 = t110 * t105;
	t108 = t106 * t107;
	t140 = t106 * t110;
	t119 = qJD(3) * (t104 * t108 + t140);
	t154 = (t96 * t119 + t127 * t143) / t94 ^ 2;
	t153 = t81 * t76;
	t152 = t96 / t102 ^ 2;
	t149 = t102 * t82;
	t147 = t87 * t110;
	t125 = 0.1e1 + t141;
	t122 = t125 * t92;
	t85 = t101 * t122;
	t146 = qJD(3) * t85;
	t139 = t107 * t110;
	t136 = 0.2e1 * t156;
	t90 = t107 * t152 + 0.1e1;
	t135 = 0.2e1 * (t108 * qJD(3) * t110 * t152 + t107 * t161) / t90 ^ 2;
	t134 = t81 * t157;
	t132 = t105 * t106 * t92;
	t130 = 0.1e1 + t152;
	t129 = t110 * t158;
	t126 = 0.2e1 * t82 * t157;
	t124 = -0.2e1 * t106 * t154;
	t123 = t101 * t132;
	t120 = t130 * t139;
	t88 = 0.1e1 / t90;
	t74 = (-t87 * t123 + t86 * t129) * t102;
	t73 = -t92 * t140 * t145 + (qJD(3) * t122 + t110 * t124) * t102;
	t72 = t85 * t148 + t147 + (-t85 * t147 - t148) * t101;
	t71 = t122 * t143 + 0.2e1 * (t119 * t92 - t125 * t154) * t101;
	t69 = -t88 * qJD(3) * t120 + (t130 * t135 - 0.2e1 * t88 * t161) * t106;
	t67 = (-t137 * t153 + (0.2e1 * t134 + (t103 * t74 + t70) * t155) * t110) * t101 + (t74 * t126 * t110 + (-t74 * t128 + (t74 * t136 + ((-t75 * t123 - t158 * t137 + t154 * t159) * t86 + (-t75 * t129 + (t105 * t124 + (t104 * t107 + t159) * t92 * qJD(3)) * t101) * t87) * t149) * t110 + (-t81 + (-(t96 - t97) * t87 * t132 + t158 * t131) * t82) * t142) * t76) * t102;
	t1 = [t73, t73, t71, 0; t67, t67, (-t145 * t153 + (-0.2e1 * t134 + (-qJD(3) * t72 - t70) * t155) * t102) * t111 + (t72 * t102 * t126 + (-t102 * qJD(3) * t81 + (t102 * t136 + t133) * t72 + (-((t71 - t143) * t86 + (t75 * t85 + qJD(3) + (-t75 - t146) * t101) * t87) * t111 + (-(-t101 * t71 - t143 * t85) * t87 - (t151 * t85 - t146 - t162) * t86) * t110) * t149) * t76) * t110, 0; t69, t69, t135 * t139 * t150 + (-t103 * t120 + (-0.2e1 * t105 * t108 - t106) * t98 * t138) * t88, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,4);
end