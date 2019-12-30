% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRP5
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
%   Wie in S5RRRRP5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:30
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRRP5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP5_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP5_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP5_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:30:30
	% EndTime: 2019-12-29 20:30:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:30:25
	% EndTime: 2019-12-29 20:30:25
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:30:25
	% EndTime: 2019-12-29 20:30:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:30:25
	% EndTime: 2019-12-29 20:30:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:30:25
	% EndTime: 2019-12-29 20:30:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:30:26
	% EndTime: 2019-12-29 20:30:27
	% DurationCPUTime: 1.52s
	% Computational Cost: add. (6770->73), mult. (3885->163), div. (1002->13), fcn. (4527->7), ass. (0->75)
	t124 = sin(qJ(1));
	t119 = t124 ^ 2;
	t117 = qJ(2) + qJ(3) + qJ(4);
	t114 = sin(t117);
	t110 = t114 ^ 2;
	t115 = cos(t117);
	t112 = 0.1e1 / t115 ^ 2;
	t160 = t110 * t112;
	t105 = t119 * t160 + 0.1e1;
	t109 = t114 * t110;
	t111 = 0.1e1 / t115;
	t113 = t111 * t112;
	t116 = qJD(2) + qJD(3) + qJD(4);
	t158 = t111 * t114;
	t133 = t116 * (t109 * t113 + t158);
	t125 = cos(qJ(1));
	t150 = qJD(1) * t125;
	t143 = t124 * t150;
	t165 = 0.1e1 / t105 ^ 2 * (t119 * t133 + t143 * t160);
	t175 = -0.2e1 * t165;
	t103 = 0.1e1 / t105;
	t140 = 0.1e1 + t160;
	t173 = t124 * t140;
	t98 = t103 * t173;
	t174 = t124 * t98 - 0.1e1;
	t121 = 0.1e1 / t125;
	t153 = t121 * t124;
	t120 = t125 ^ 2;
	t172 = qJD(1) * (t119 / t120 + 0.1e1) * t153;
	t152 = t124 * t114;
	t102 = atan2(-t152, -t115);
	t101 = cos(t102);
	t100 = sin(t102);
	t145 = t100 * t152;
	t97 = -t101 * t115 - t145;
	t94 = 0.1e1 / t97;
	t95 = 0.1e1 / t97 ^ 2;
	t171 = t103 - 0.1e1;
	t156 = t115 * t116;
	t146 = t95 * t156;
	t142 = t114 * t150;
	t155 = t116 * t124;
	t161 = t101 * t114;
	t144 = t112 * t155;
	t89 = (-(-t115 * t155 - t142) * t111 + t110 * t144) * t103;
	t84 = (-t124 * t89 + t116) * t161 + (-t142 + (t89 - t155) * t115) * t100;
	t169 = t84 * t94 * t95;
	t92 = t120 * t110 * t95 + 0.1e1;
	t170 = (t120 * t114 * t146 + (-t120 * t169 - t95 * t143) * t110) / t92 ^ 2;
	t90 = 0.1e1 / t92;
	t168 = t90 * t95;
	t167 = t94 * t90;
	t164 = t116 * t98;
	t162 = t125 * t95;
	t159 = t110 * t124;
	t157 = t114 * t125;
	t154 = t119 / t125 ^ 2;
	t151 = qJD(1) * t124;
	t149 = 0.2e1 * t169;
	t108 = t112 * t154 + 0.1e1;
	t148 = 0.2e1 / t108 ^ 2 * (t113 * t116 * t114 * t154 + t112 * t172);
	t147 = t94 * t170;
	t141 = 0.2e1 * t95 * t170;
	t139 = 0.1e1 + t154;
	t138 = t111 * t175;
	t137 = t101 * t103 * t110 * t111;
	t136 = t140 * t125;
	t134 = t139 * t114 * t112;
	t106 = 0.1e1 / t108;
	t88 = (t171 * t114 * t100 - t124 * t137) * t125;
	t87 = t121 * t112 * t148 * t152 + ((-0.2e1 * t110 * t113 - t111) * t116 * t153 - qJD(1) * t134) * t106;
	t86 = -t174 * t161 + (-t124 + t98) * t115 * t100;
	t85 = t173 * t175 + (qJD(1) * t136 + 0.2e1 * t124 * t133) * t103;
	t82 = (-t151 * t167 + (-0.2e1 * t147 + (-t116 * t86 - t84) * t168) * t125) * t115 + (t86 * t125 * t141 + (-t125 * t116 * t94 - ((-t124 * t85 - t150 * t98) * t101 + (t174 * t89 + t155 - t164) * t100) * t95 * t157 + (t125 * t149 + t95 * t151) * t86 - ((t85 - t150) * t100 + (t89 * t98 + t116 + (-t89 - t164) * t124) * t101) * t115 * t162) * t90) * t114;
	t1 = [t138 * t157 + (t116 * t136 - t151 * t158) * t103, t85, t85, t85, 0; (-t156 * t167 + (0.2e1 * t147 + (qJD(1) * t88 + t84) * t168) * t114) * t124 + (t88 * t141 * t114 + (-t88 * t146 + (t88 * t149 + ((0.2e1 * t114 * t165 + t156 + (-t111 * t89 * t159 - t156) * t103) * t100 + (t138 * t159 + t114 * t89 + (t109 * t144 - (t89 - 0.2e1 * t155) * t114) * t103) * t101) * t162) * t114 + (-t94 + (-(t119 - t120) * t137 + t171 * t145) * t95) * t114 * qJD(1)) * t90) * t125, t82, t82, t82, 0; t139 * t111 * t148 + (-0.2e1 * t111 * t172 - t116 * t134) * t106, t87, t87, t87, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end