% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6RPPPRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:25
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPPRR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:52
	% EndTime: 2019-10-09 23:25:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:52
	% EndTime: 2019-10-09 23:25:52
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:52
	% EndTime: 2019-10-09 23:25:52
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:52
	% EndTime: 2019-10-09 23:25:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:52
	% EndTime: 2019-10-09 23:25:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:52
	% EndTime: 2019-10-09 23:25:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:53
	% EndTime: 2019-10-09 23:25:54
	% DurationCPUTime: 0.96s
	% Computational Cost: add. (1598->92), mult. (2519->207), div. (480->12), fcn. (2968->9), ass. (0->93)
	t124 = qJ(1) + pkin(9);
	t122 = sin(t124);
	t131 = sin(qJ(5));
	t126 = 0.1e1 / t131 ^ 2;
	t133 = cos(qJ(5));
	t129 = t133 ^ 2;
	t176 = t126 * t129;
	t152 = 0.1e1 + t176;
	t194 = t122 * t152;
	t130 = sin(qJ(6));
	t132 = cos(qJ(6));
	t149 = qJD(6) * t131 + qJD(1);
	t168 = qJD(5) * t133;
	t193 = t149 * t130 - t132 * t168;
	t179 = t122 * t133;
	t115 = atan2(t179, t131);
	t112 = cos(t115);
	t111 = sin(t115);
	t162 = t111 * t179;
	t102 = t112 * t131 + t162;
	t99 = 0.1e1 / t102;
	t123 = cos(t124);
	t174 = t131 * t132;
	t160 = t123 * t174;
	t180 = t122 * t130;
	t110 = t160 - t180;
	t104 = 0.1e1 / t110;
	t125 = 0.1e1 / t131;
	t100 = 0.1e1 / t102 ^ 2;
	t105 = 0.1e1 / t110 ^ 2;
	t120 = t122 ^ 2;
	t116 = t120 * t176 + 0.1e1;
	t113 = 0.1e1 / t116;
	t192 = t113 - 0.1e1;
	t175 = t130 * t131;
	t109 = t122 * t132 + t123 * t175;
	t103 = t109 ^ 2;
	t184 = t105 * t109;
	t148 = qJD(1) * t131 + qJD(6);
	t143 = t148 * t132;
	t90 = -t122 * t143 - t193 * t123;
	t188 = t104 * t105 * t90;
	t155 = t130 * t168;
	t172 = qJD(1) * t123;
	t89 = -qJD(6) * t160 - t123 * t155 - t132 * t172 + t148 * t180;
	t98 = t103 * t105 + 0.1e1;
	t191 = (-t103 * t188 - t89 * t184) / t98 ^ 2;
	t121 = t123 ^ 2;
	t182 = t121 * t129;
	t93 = t100 * t182 + 0.1e1;
	t91 = 0.1e1 / t93;
	t190 = t100 * t91;
	t171 = qJD(1) * t133;
	t157 = t123 * t171;
	t170 = qJD(5) * t122;
	t183 = t112 * t133;
	t156 = t126 * t170;
	t169 = qJD(5) * t131;
	t88 = ((-t122 * t169 + t157) * t125 - t129 * t156) * t113;
	t83 = (t122 * t88 + qJD(5)) * t183 + (t157 + (-t88 - t170) * t131) * t111;
	t189 = t99 * t100 * t83;
	t128 = t133 * t129;
	t144 = (t126 * t128 + t133) * t125;
	t181 = t122 * t129;
	t146 = t172 * t181;
	t187 = 0.1e1 / t116 ^ 2 * (-t144 * t120 * qJD(5) + t126 * t146);
	t97 = t113 * t194;
	t186 = -t122 + t97;
	t178 = t123 * t130;
	t177 = t123 * t133;
	t173 = qJD(1) * t122;
	t167 = -0.2e1 * (-t182 * t189 + (-t121 * t131 * t168 - t146) * t100) / t93 ^ 2;
	t166 = 0.2e1 * t191;
	t165 = -0.2e1 * t189;
	t164 = 0.2e1 * t187;
	t163 = t91 * t169;
	t161 = t113 * t125 * t129;
	t159 = -t122 * t97 + 0.1e1;
	t158 = t122 * t171;
	t153 = t99 * t167;
	t151 = 0.2e1 * t109 * t188;
	t150 = -0.2e1 * t125 * t187;
	t147 = t122 * t161;
	t145 = t152 * t123;
	t142 = -t104 * t130 + t132 * t184;
	t95 = 0.1e1 / t98;
	t141 = t142 * t95;
	t108 = -t122 * t174 - t178;
	t107 = -t122 * t175 + t123 * t132;
	t87 = (-t192 * t133 * t111 + t112 * t147) * t123;
	t86 = t186 * t131 * t111 + t159 * t183;
	t84 = t164 * t194 + (-qJD(1) * t145 + 0.2e1 * t144 * t170) * t113;
	t1 = [t150 * t177 + (-qJD(5) * t145 - t125 * t158) * t113, 0, 0, 0, t84, 0; (-t99 * t163 + (t153 + (-qJD(1) * t87 - t83) * t190) * t133) * t122 + (t87 * t133 * t91 * t165 + (-t87 * t163 + (t87 * t167 + ((t133 * t164 - t88 * t147 + t192 * t169) * t111 + (t150 * t181 + t133 * t88 + (-t128 * t156 + (-t88 - 0.2e1 * t170) * t133) * t113) * t112) * t91 * t123) * t133) * t100 + (t99 + ((-t120 + t121) * t112 * t161 + t192 * t162) * t100) * t91 * t171) * t123, 0, 0, 0, (-t99 * t91 * t173 + (t153 + (-qJD(5) * t86 - t83) * t190) * t123) * t131 + (((qJD(5) * t99 + t86 * t165) * t123 + (-t86 * t173 + ((t122 * t84 - t97 * t172) * t112 + (t186 * qJD(5) - t159 * t88) * t111) * t177) * t100) * t91 + (t86 * t167 + ((-t84 - t172) * t111 + (t88 * t97 - qJD(5) + (qJD(5) * t97 - t88) * t122) * t112) * t91 * t131) * t100 * t123) * t133, 0; (-t104 * t107 + t108 * t184) * t166 + (t108 * t151 + (-t107 * t90 + t108 * t89 + (-t193 * t122 + t123 * t143) * t109) * t105 + (-t148 * t178 + (-t149 * t132 - t155) * t122) * t104) * t95, 0, 0, 0, t141 * t158 + (t141 * t169 + (t142 * t166 + ((qJD(6) * t104 + t151) * t132 + (t132 * t89 + (qJD(6) * t109 - t90) * t130) * t105) * t95) * t133) * t123, -0.2e1 * t191 + 0.2e1 * (-t105 * t89 * t95 + (-t105 * t191 - t95 * t188) * t109) * t109;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end