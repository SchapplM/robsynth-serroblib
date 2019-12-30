% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRRP9
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
%   Wie in S5RPRRP9_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPRRP9_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP9_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:27:17
	% EndTime: 2019-12-29 17:27:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:27:12
	% EndTime: 2019-12-29 17:27:12
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:27:18
	% EndTime: 2019-12-29 17:27:18
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:27:12
	% EndTime: 2019-12-29 17:27:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:27:12
	% EndTime: 2019-12-29 17:27:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:27:13
	% EndTime: 2019-12-29 17:27:14
	% DurationCPUTime: 1.36s
	% Computational Cost: add. (4829->71), mult. (2860->161), div. (736->13), fcn. (3352->7), ass. (0->75)
	t111 = pkin(8) + qJ(3) + qJ(4);
	t110 = cos(t111);
	t106 = 0.1e1 / t110;
	t170 = -0.2e1 * t106;
	t119 = sin(qJ(1));
	t109 = sin(t111);
	t105 = t109 ^ 2;
	t107 = 0.1e1 / t110 ^ 2;
	t153 = t105 * t107;
	t135 = 0.1e1 + t153;
	t114 = t119 ^ 2;
	t100 = t114 * t153 + 0.1e1;
	t98 = 0.1e1 / t100;
	t131 = t135 * t98;
	t93 = t119 * t131;
	t169 = t119 * t93 - 0.1e1;
	t120 = cos(qJ(1));
	t116 = 0.1e1 / t120;
	t149 = t116 * t119;
	t115 = t120 ^ 2;
	t168 = qJD(1) * (t114 / t115 + 0.1e1) * t149;
	t148 = t119 * t109;
	t97 = atan2(-t148, -t110);
	t95 = sin(t97);
	t139 = t95 * t148;
	t96 = cos(t97);
	t92 = -t96 * t110 - t139;
	t89 = 0.1e1 / t92;
	t90 = 0.1e1 / t92 ^ 2;
	t167 = 0.2e1 * t109;
	t166 = t98 - 0.1e1;
	t145 = qJD(1) * t120;
	t137 = t119 * t145;
	t112 = qJD(3) + qJD(4);
	t152 = t110 * t112;
	t140 = t90 * t152;
	t151 = t112 * t119;
	t157 = t110 * t95;
	t84 = (-(-t109 * t145 - t110 * t151) * t106 + t151 * t153) * t98;
	t79 = (t84 - t151) * t157 + (-t95 * t145 + (-t119 * t84 + t112) * t96) * t109;
	t164 = t79 * t89 * t90;
	t159 = t105 * t90;
	t87 = t115 * t159 + 0.1e1;
	t165 = (t115 * t109 * t140 + (-t115 * t164 - t90 * t137) * t105) / t87 ^ 2;
	t85 = 0.1e1 / t87;
	t162 = t85 * t90;
	t104 = t109 * t105;
	t108 = t106 * t107;
	t128 = t112 * (t104 * t108 + t106 * t109);
	t161 = (t114 * t128 + t137 * t153) / t100 ^ 2;
	t160 = t89 * t85;
	t158 = t106 * t98;
	t156 = t112 * t93;
	t154 = t120 * t90;
	t150 = t114 / t120 ^ 2;
	t147 = qJD(1) * t109;
	t146 = qJD(1) * t119;
	t144 = 0.2e1 * t164;
	t103 = t107 * t150 + 0.1e1;
	t143 = 0.2e1 / t103 ^ 2 * (t108 * t112 * t109 * t150 + t107 * t168);
	t142 = t89 * t165;
	t141 = t119 * t158;
	t138 = t166 * t109;
	t136 = 0.2e1 * t90 * t165;
	t134 = 0.1e1 + t150;
	t133 = t161 * t170;
	t132 = t105 * t141;
	t129 = t134 * t109 * t107;
	t101 = 0.1e1 / t103;
	t83 = (-t96 * t132 + t95 * t138) * t120;
	t82 = t116 * t107 * t143 * t148 + ((-0.2e1 * t105 * t108 - t106) * t112 * t149 - qJD(1) * t129) * t101;
	t81 = (-t119 + t93) * t157 - t169 * t96 * t109;
	t80 = t131 * t145 + 0.2e1 * (t128 * t98 - t135 * t161) * t119;
	t77 = (-t146 * t160 + (-0.2e1 * t142 + (-t112 * t81 - t79) * t162) * t120) * t110 + (t81 * t120 * t136 + (-t120 * t112 * t89 + (t120 * t144 + t90 * t146) * t81 + (-((-t119 * t80 - t145 * t93) * t96 + (t169 * t84 + t151 - t156) * t95) * t109 - ((t80 - t145) * t95 + (t84 * t93 + t112 + (-t84 - t156) * t119) * t96) * t110) * t154) * t85) * t109;
	t1 = [-t141 * t147 + (t109 * t133 + t112 * t131) * t120, 0, t80, t80, 0; (-t152 * t160 + (0.2e1 * t142 + (qJD(1) * t83 + t79) * t162) * t109) * t119 + (t83 * t136 * t109 + (-t83 * t140 + (t83 * t144 + ((-t84 * t132 - t166 * t152 + t161 * t167) * t95 + (-t84 * t138 + (t105 * t133 + (t104 * t107 + t167) * t98 * t112) * t119) * t96) * t154) * t109 + (-t89 + t166 * t90 * t139 - (t114 - t115) * t96 * t158 * t159) * t147) * t85) * t120, 0, t77, t77, 0; t134 * t106 * t143 + (-t112 * t129 + t168 * t170) * t101, 0, t82, t82, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end