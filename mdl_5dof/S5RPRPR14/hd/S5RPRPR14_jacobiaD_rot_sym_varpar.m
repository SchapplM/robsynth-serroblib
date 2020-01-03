% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRPR14
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
%   Wie in S5RPRPR14_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPRPR14_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR14_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:35:47
	% EndTime: 2019-12-31 18:35:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:35:47
	% EndTime: 2019-12-31 18:35:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:35:47
	% EndTime: 2019-12-31 18:35:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:35:47
	% EndTime: 2019-12-31 18:35:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:35:47
	% EndTime: 2019-12-31 18:35:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:35:47
	% EndTime: 2019-12-31 18:35:48
	% DurationCPUTime: 0.66s
	% Computational Cost: add. (2081->92), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->95)
	t139 = cos(qJ(1));
	t203 = 0.2e1 * t139;
	t133 = qJ(3) + pkin(8);
	t131 = sin(t133);
	t132 = cos(t133);
	t181 = t139 * t132;
	t121 = atan2(-t181, t131);
	t119 = sin(t121);
	t120 = cos(t121);
	t105 = -t119 * t181 + t120 * t131;
	t102 = 0.1e1 / t105;
	t136 = sin(qJ(5));
	t180 = t139 * t136;
	t137 = sin(qJ(1));
	t138 = cos(qJ(5));
	t182 = t137 * t138;
	t116 = t131 * t182 + t180;
	t112 = 0.1e1 / t116;
	t126 = 0.1e1 / t131;
	t103 = 0.1e1 / t105 ^ 2;
	t113 = 0.1e1 / t116 ^ 2;
	t127 = 0.1e1 / t131 ^ 2;
	t135 = t139 ^ 2;
	t130 = t132 ^ 2;
	t186 = t127 * t130;
	t124 = t135 * t186 + 0.1e1;
	t122 = 0.1e1 / t124;
	t202 = t122 - 0.1e1;
	t134 = t137 ^ 2;
	t177 = qJD(1) * t139;
	t155 = t130 * t137 * t177;
	t175 = qJD(3) * t132;
	t185 = t130 * t134;
	t174 = qJD(3) * t139;
	t165 = t127 * t174;
	t178 = qJD(1) * t137;
	t166 = t132 * t178;
	t96 = ((t131 * t174 + t166) * t126 + t130 * t165) * t122;
	t160 = -t96 + t174;
	t161 = -t139 * t96 + qJD(3);
	t189 = t120 * t132;
	t91 = t161 * t189 + (t160 * t131 + t166) * t119;
	t199 = t102 * t103 * t91;
	t99 = t103 * t185 + 0.1e1;
	t201 = (-t185 * t199 + (-t131 * t134 * t175 + t155) * t103) / t99 ^ 2;
	t97 = 0.1e1 / t99;
	t200 = t103 * t97;
	t158 = qJD(1) * t131 + qJD(5);
	t151 = t158 * t139;
	t159 = qJD(5) * t131 + qJD(1);
	t152 = t159 * t138;
	t100 = t137 * t152 + (t137 * t175 + t151) * t136;
	t179 = t139 * t138;
	t183 = t137 * t136;
	t115 = t131 * t183 - t179;
	t111 = t115 ^ 2;
	t110 = t111 * t113 + 0.1e1;
	t192 = t113 * t115;
	t153 = t159 * t136;
	t101 = t138 * t151 + (t138 * t175 - t153) * t137;
	t196 = t101 * t112 * t113;
	t198 = 0.1e1 / t110 ^ 2 * (t100 * t192 - t111 * t196);
	t129 = t132 * t130;
	t187 = t126 * t132;
	t149 = qJD(3) * (-t126 * t127 * t129 - t187);
	t194 = (-t127 * t155 + t135 * t149) / t124 ^ 2;
	t193 = t112 * t136;
	t191 = t115 * t138;
	t190 = t119 * t131;
	t188 = t126 * t130;
	t176 = qJD(3) * t131;
	t173 = -0.2e1 * t201;
	t172 = -0.2e1 * t199;
	t171 = 0.2e1 * t198;
	t170 = t102 * t201;
	t169 = t97 * t176;
	t168 = t132 * t194;
	t167 = t122 * t188;
	t164 = 0.1e1 + t186;
	t163 = t194 * t203;
	t162 = 0.2e1 * t115 * t196;
	t157 = t139 * t167;
	t156 = t202 * t132 * t119;
	t154 = t164 * t137;
	t150 = t113 * t191 - t193;
	t148 = t150 * t137;
	t147 = t132 * t174 - t158 * t137;
	t118 = t131 * t179 - t183;
	t117 = t131 * t180 + t182;
	t108 = 0.1e1 / t110;
	t107 = t164 * t139 * t122;
	t95 = (-t120 * t157 - t156) * t137;
	t93 = t139 * t190 + t189 + (-t120 * t181 - t190) * t107;
	t92 = -t164 * t163 + (-qJD(1) * t154 + t149 * t203) * t122;
	t1 = [-0.2e1 * t137 * t126 * t168 + (-qJD(3) * t154 + t177 * t187) * t122, 0, t92, 0, 0; (t102 * t169 + (0.2e1 * t170 + (qJD(1) * t95 + t91) * t200) * t132) * t139 + ((-t95 * t169 + (t95 * t173 + ((t96 * t157 + t202 * t176 + 0.2e1 * t168) * t119 + (t163 * t188 + t132 * t96 + (t129 * t165 + (-t96 + 0.2e1 * t174) * t132) * t122) * t120) * t97 * t137) * t132) * t103 + (t95 * t172 + (t102 + ((t134 - t135) * t120 * t167 - t139 * t156) * t103) * qJD(1)) * t132 * t97) * t137, 0, (t102 * t97 * t177 + (-0.2e1 * t170 + (-qJD(3) * t93 - t91) * t200) * t137) * t131 + (((qJD(3) * t102 + t93 * t172) * t137 + (t93 * t177 + ((t107 * t178 - t139 * t92) * t120 + ((t107 * t139 - 0.1e1) * t96 + (-t107 + t139) * qJD(3)) * t119) * t132 * t137) * t103) * t97 + (t93 * t173 + ((-t92 - t178) * t119 + (t160 * t107 - t161) * t120) * t97 * t131) * t103 * t137) * t132, 0, 0; (-t112 * t117 + t118 * t192) * t171 + (t118 * t162 + t139 * t112 * t152 + t147 * t193 + (t139 * t115 * t153 - t118 * t100 - t117 * t101 - t147 * t191) * t113) * t108, 0, t132 * t148 * t171 + (t148 * t176 + (-t150 * t177 + ((qJD(5) * t112 + t162) * t138 + (-t100 * t138 + (qJD(5) * t115 - t101) * t136) * t113) * t137) * t132) * t108, 0, -0.2e1 * t198 + 0.2e1 * (t100 * t113 * t108 + (-t108 * t196 - t113 * t198) * t115) * t115;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end