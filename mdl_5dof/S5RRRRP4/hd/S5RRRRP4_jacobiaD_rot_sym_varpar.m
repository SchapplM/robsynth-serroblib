% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRP4
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RRRRP4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRRP4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:27:40
	% EndTime: 2019-12-29 20:27:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:27:40
	% EndTime: 2019-12-29 20:27:40
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:27:40
	% EndTime: 2019-12-29 20:27:40
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (74->0), mult. (74->0), div. (30->0), fcn. (44->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:27:45
	% EndTime: 2019-12-29 20:27:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:27:40
	% EndTime: 2019-12-29 20:27:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:27:45
	% EndTime: 2019-12-29 20:27:47
	% DurationCPUTime: 1.43s
	% Computational Cost: add. (5514->75), mult. (3670->165), div. (940->13), fcn. (4354->7), ass. (0->81)
	t139 = qJ(1) + qJ(2);
	t133 = sin(t139);
	t138 = qJ(3) + qJ(4);
	t132 = sin(t138);
	t122 = t132 ^ 2;
	t134 = cos(t138);
	t126 = 0.1e1 / t134 ^ 2;
	t179 = t122 * t126;
	t157 = 0.1e1 + t179;
	t187 = t133 * t157;
	t173 = t133 * t132;
	t114 = atan2(-t173, -t134);
	t112 = sin(t114);
	t113 = cos(t114);
	t110 = -t112 * t173 - t113 * t134;
	t107 = 0.1e1 / t110;
	t125 = 0.1e1 / t134;
	t135 = cos(t139);
	t129 = 0.1e1 / t135;
	t108 = 0.1e1 / t110 ^ 2;
	t124 = t133 ^ 2;
	t119 = t124 * t179 + 0.1e1;
	t115 = 0.1e1 / t119;
	t137 = qJD(1) + qJD(2);
	t169 = t135 * t137;
	t158 = t132 * t169;
	t136 = qJD(3) + qJD(4);
	t172 = t133 * t136;
	t160 = t126 * t172;
	t170 = t134 * t136;
	t101 = (-(-t133 * t170 - t158) * t125 + t122 * t160) * t115;
	t153 = t101 * t133 - t136;
	t154 = t101 - t172;
	t180 = t113 * t132;
	t94 = -t153 * t180 + (t154 * t134 - t158) * t112;
	t186 = t107 * t108 * t94;
	t177 = t122 * t133;
	t161 = t125 * t177;
	t152 = t113 * t161;
	t100 = (-t115 * t152 + (t115 - 0.1e1) * t132 * t112) * t135;
	t185 = t100 * t108;
	t121 = t132 * t122;
	t127 = t125 * t126;
	t147 = (t121 * t127 + t125 * t132) * t136;
	t151 = t169 * t177;
	t184 = (t124 * t147 + t126 * t151) / t119 ^ 2;
	t183 = t108 * t132;
	t182 = t108 * t135;
	t181 = t112 * t134;
	t128 = t135 ^ 2;
	t178 = t122 * t128;
	t176 = t124 / t135 ^ 2;
	t175 = t125 * t137;
	t174 = t132 * t135;
	t171 = t133 * t137;
	t168 = t136 * t107;
	t167 = t136 * t132;
	t166 = t137 * t107;
	t104 = t108 * t178 + 0.1e1;
	t165 = 0.2e1 / t104 ^ 2 * (-t178 * t186 + (t128 * t134 * t167 - t151) * t108);
	t164 = 0.2e1 * t186;
	t163 = -0.2e1 * t184;
	t120 = t126 * t176 + 0.1e1;
	t149 = (t124 / t128 + 0.1e1) * t129 * t133;
	t162 = 0.2e1 * (t149 * t137 * t126 + t127 * t167 * t176) / t120 ^ 2;
	t159 = t132 * t171;
	t156 = 0.1e1 + t176;
	t155 = t132 * t165;
	t150 = t135 * t157;
	t148 = t126 * t132 * t156;
	t117 = 0.1e1 / t120;
	t111 = t115 * t187;
	t102 = 0.1e1 / t104;
	t99 = t125 * t163 * t174 + (-t125 * t159 + t136 * t150) * t115;
	t98 = -t133 * t181 + t180 + (-t113 * t173 + t181) * t111;
	t97 = t129 * t126 * t162 * t173 + ((-0.2e1 * t122 * t127 - t125) * t129 * t172 - t137 * t148) * t117;
	t96 = t163 * t187 + (0.2e1 * t133 * t147 + t137 * t150) * t115;
	t95 = t156 * t125 * t162 + (-t136 * t148 - 0.2e1 * t149 * t175) * t117;
	t92 = (t107 * t155 + (-t134 * t168 + (t100 * t137 + t94) * t183) * t102) * t133 + (t155 * t185 + (-t170 * t185 + (t100 * t164 - t166 + (-(-t112 * t170 + 0.2e1 * t152 * t184) * t135 + (-t112 * t171 - (-t101 * t113 + t112 * t163) * t135) * t132 + ((-(t101 * t161 + t170) * t135 + t159) * t112 + ((t121 * t160 - (t101 - 0.2e1 * t172) * t132) * t135 + (-t124 + t128) * t122 * t175) * t113) * t115) * t108) * t132) * t102) * t135;
	t91 = (-t107 * t134 + t98 * t183) * t135 * t165 + ((-t133 * t166 + (-t136 * t98 - t94) * t182) * t134 + ((t98 * t164 - t168) * t135 + (t98 * t171 + (-(-t111 * t169 - t133 * t96) * t113 - (t153 * t111 - t154) * t112) * t174) * t108 - ((t96 - t169) * t112 + ((-t111 * t133 + 0.1e1) * t136 + (t111 - t133) * t101) * t113) * t134 * t182) * t132) * t102;
	t1 = [t99, t99, t96, t96, 0; t92, t92, t91, t91, 0; t95, t95, t97, t97, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end