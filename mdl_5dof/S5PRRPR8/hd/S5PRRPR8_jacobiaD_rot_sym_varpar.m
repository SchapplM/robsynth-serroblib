% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRPR8
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
%   Wie in S5PRRPR8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PRRPR8_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR8_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:40:16
	% EndTime: 2019-12-29 15:40:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:40:21
	% EndTime: 2019-12-29 15:40:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:40:16
	% EndTime: 2019-12-29 15:40:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:40:16
	% EndTime: 2019-12-29 15:40:16
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:40:21
	% EndTime: 2019-12-29 15:40:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:40:22
	% EndTime: 2019-12-29 15:40:23
	% DurationCPUTime: 0.82s
	% Computational Cost: add. (3038->45), mult. (2152->102), div. (440->12), fcn. (2434->9), ass. (0->60)
	t141 = qJ(2) + qJ(3) + pkin(9);
	t139 = sin(t141);
	t134 = t139 ^ 2;
	t140 = cos(t141);
	t170 = t134 / t140 ^ 2;
	t160 = 0.1e1 + t170;
	t144 = qJD(2) + qJD(3);
	t145 = sin(pkin(8));
	t180 = t144 * t145;
	t147 = sin(qJ(5));
	t148 = cos(qJ(5));
	t146 = cos(pkin(8));
	t168 = t140 * t146;
	t126 = t145 * t147 + t148 * t168;
	t167 = t145 * t139;
	t127 = atan2(-t167, -t140);
	t119 = sin(t127);
	t120 = cos(t127);
	t111 = -t119 * t167 - t120 * t140;
	t108 = 0.1e1 / t111;
	t122 = 0.1e1 / t126;
	t109 = 0.1e1 / t111 ^ 2;
	t123 = 0.1e1 / t126 ^ 2;
	t142 = t145 ^ 2;
	t130 = t142 * t170 + 0.1e1;
	t128 = 0.1e1 / t130;
	t112 = t160 * t145 * t128;
	t172 = t119 * t140;
	t158 = t120 * t139 - t145 * t172;
	t163 = t120 * t167;
	t159 = -t163 + t172;
	t101 = t159 * t112 + t158;
	t179 = 0.2e1 * t101;
	t143 = t146 ^ 2;
	t106 = t143 * t134 * t109 + 0.1e1;
	t169 = t140 * t144;
	t174 = t109 * t139;
	t107 = t112 * t144;
	t100 = t159 * t107 + t158 * t144;
	t177 = t100 * t108 * t109;
	t178 = 0.1e1 / t106 ^ 2 * (-t134 * t177 + t169 * t174) * t143;
	t125 = -t145 * t148 + t147 * t168;
	t121 = t125 ^ 2;
	t115 = t121 * t123 + 0.1e1;
	t162 = t139 * t144 * t146;
	t117 = -t126 * qJD(5) + t147 * t162;
	t171 = t123 * t125;
	t165 = qJD(5) * t125;
	t118 = -t148 * t162 - t165;
	t173 = t118 * t122 * t123;
	t176 = (-t117 * t171 - t121 * t173) / t115 ^ 2;
	t104 = 0.1e1 / t106;
	t175 = t104 * t109;
	t164 = -0.2e1 * t176;
	t157 = -t122 * t147 + t148 * t171;
	t113 = 0.1e1 / t115;
	t102 = 0.2e1 * (t128 - t160 / t130 ^ 2 * t142) * t160 / t140 * t139 * t180;
	t98 = (t157 * t139 * t164 + (t157 * t169 + ((-qJD(5) * t122 - 0.2e1 * t125 * t173) * t148 + (-t117 * t148 + (t118 - t165) * t147) * t123) * t139) * t113) * t146;
	t97 = ((-0.2e1 * t108 * t178 + (-t101 * t144 - t100) * t175) * t140 + (t109 * t178 * t179 + (t102 * t109 * t163 + t177 * t179 - t144 * t108 - (t180 - t107 + (t107 * t145 - t144) * t112) * t119 * t174) * t104 - (t102 * t119 + ((-t112 * t145 + 0.1e1) * t144 + (t112 - t145) * t107) * t120) * t140 * t175) * t139) * t146;
	t1 = [0, t102, t102, 0, 0; 0, t97, t97, 0, 0; 0, t98, t98, 0, t164 + 0.2e1 * (-t113 * t117 * t123 + (-t113 * t173 - t123 * t176) * t125) * t125;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end