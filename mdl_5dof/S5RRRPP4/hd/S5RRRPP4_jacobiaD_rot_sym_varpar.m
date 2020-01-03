% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRPP4
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
%   Wie in S5RRRPP4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRPP4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:56:45
	% EndTime: 2019-12-31 20:56:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:56:45
	% EndTime: 2019-12-31 20:56:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:56:45
	% EndTime: 2019-12-31 20:56:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:56:45
	% EndTime: 2019-12-31 20:56:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:56:45
	% EndTime: 2019-12-31 20:56:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:56:45
	% EndTime: 2019-12-31 20:56:46
	% DurationCPUTime: 0.60s
	% Computational Cost: add. (4829->74), mult. (2860->166), div. (736->13), fcn. (3352->7), ass. (0->75)
	t123 = sin(qJ(1));
	t118 = t123 ^ 2;
	t115 = qJ(2) + qJ(3) + pkin(8);
	t113 = sin(t115);
	t109 = t113 ^ 2;
	t114 = cos(t115);
	t111 = 0.1e1 / t114 ^ 2;
	t158 = t109 * t111;
	t104 = t118 * t158 + 0.1e1;
	t108 = t113 * t109;
	t110 = 0.1e1 / t114;
	t112 = t110 * t111;
	t116 = qJD(2) + qJD(3);
	t156 = t110 * t113;
	t132 = t116 * (t108 * t112 + t156);
	t124 = cos(qJ(1));
	t148 = qJD(1) * t124;
	t140 = t123 * t148;
	t165 = 0.1e1 / t104 ^ 2 * (t118 * t132 + t140 * t158);
	t174 = -0.2e1 * t165;
	t102 = 0.1e1 / t104;
	t138 = 0.1e1 + t158;
	t172 = t123 * t138;
	t97 = t102 * t172;
	t173 = t123 * t97 - 0.1e1;
	t120 = 0.1e1 / t124;
	t151 = t120 * t123;
	t119 = t124 ^ 2;
	t171 = qJD(1) * (t118 / t119 + 0.1e1) * t151;
	t150 = t123 * t113;
	t101 = atan2(-t150, -t114);
	t100 = cos(t101);
	t99 = sin(t101);
	t143 = t99 * t150;
	t96 = -t100 * t114 - t143;
	t93 = 0.1e1 / t96;
	t94 = 0.1e1 / t96 ^ 2;
	t154 = t114 * t116;
	t144 = t94 * t154;
	t153 = t116 * t123;
	t162 = t114 * t99;
	t141 = t111 * t153;
	t88 = (-(-t113 * t148 - t114 * t153) * t110 + t109 * t141) * t102;
	t83 = (t88 - t153) * t162 + (-t99 * t148 + (-t123 * t88 + t116) * t100) * t113;
	t169 = t83 * t93 * t94;
	t164 = t109 * t94;
	t91 = t119 * t164 + 0.1e1;
	t170 = (t119 * t113 * t144 + (-t119 * t169 - t94 * t140) * t109) / t91 ^ 2;
	t89 = 0.1e1 / t91;
	t167 = t89 * t94;
	t166 = t93 * t89;
	t163 = t113 * t99;
	t161 = t116 * t97;
	t159 = t124 * t94;
	t157 = t109 * t123;
	t155 = t113 * t124;
	t152 = t118 / t124 ^ 2;
	t149 = qJD(1) * t123;
	t147 = 0.2e1 * t169;
	t107 = t111 * t152 + 0.1e1;
	t146 = 0.2e1 / t107 ^ 2 * (t112 * t116 * t113 * t152 + t111 * t171);
	t145 = t93 * t170;
	t142 = t110 * t157;
	t139 = 0.2e1 * t94 * t170;
	t137 = 0.1e1 + t152;
	t136 = t110 * t174;
	t135 = t138 * t124;
	t133 = t137 * t113 * t111;
	t105 = 0.1e1 / t107;
	t87 = (-t163 + (-t100 * t142 + t163) * t102) * t124;
	t86 = t120 * t111 * t146 * t150 + ((-0.2e1 * t109 * t112 - t110) * t116 * t151 - qJD(1) * t133) * t105;
	t85 = (-t123 + t97) * t162 - t173 * t113 * t100;
	t84 = t172 * t174 + (qJD(1) * t135 + 0.2e1 * t123 * t132) * t102;
	t81 = (-t149 * t166 + (-0.2e1 * t145 + (-t116 * t85 - t83) * t167) * t124) * t114 + (t85 * t124 * t139 + (-t124 * t116 * t93 - ((-t123 * t84 - t148 * t97) * t100 + (t173 * t88 + t153 - t161) * t99) * t94 * t155 + (t124 * t147 + t94 * t149) * t85 - ((t84 - t148) * t99 + (t88 * t97 + t116 + (-t88 - t161) * t123) * t100) * t114 * t159) * t89) * t113;
	t1 = [t136 * t155 + (t116 * t135 - t149 * t156) * t102, t84, t84, 0, 0; (-t154 * t166 + (0.2e1 * t145 + (qJD(1) * t87 + t83) * t167) * t113) * t123 + (t87 * t139 * t113 + (-t87 * t144 + (t87 * t147 + ((0.2e1 * t113 * t165 + t154 + (-t88 * t142 - t154) * t102) * t99 + (t136 * t157 + t113 * t88 + (t108 * t141 - (t88 - 0.2e1 * t153) * t113) * t102) * t100) * t159) * t113 + (-t93 - (-t102 + 0.1e1) * t94 * t143 - (t118 - t119) * t110 * t102 * t100 * t164) * t113 * qJD(1)) * t89) * t124, t81, t81, 0, 0; t137 * t110 * t146 + (-0.2e1 * t110 * t171 - t116 * t133) * t105, t86, t86, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end