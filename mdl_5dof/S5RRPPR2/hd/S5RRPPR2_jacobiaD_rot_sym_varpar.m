% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPPR2
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
%   Wie in S5RRPPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRPPR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
JaD_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:06:21
	% EndTime: 2022-01-20 10:06:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:06:21
	% EndTime: 2022-01-20 10:06:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:06:22
	% EndTime: 2022-01-20 10:06:22
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (74->0), mult. (74->0), div. (30->0), fcn. (44->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:06:21
	% EndTime: 2022-01-20 10:06:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (118->0), mult. (74->0), div. (30->0), fcn. (44->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:06:22
	% EndTime: 2022-01-20 10:06:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:06:22
	% EndTime: 2022-01-20 10:06:22
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (1990->41), mult. (1512->108), div. (234->12), fcn. (1884->9), ass. (0->57)
	t124 = qJ(1) + qJ(2) + pkin(8);
	t123 = cos(t124);
	t121 = t123 ^ 2;
	t131 = sin(pkin(9));
	t126 = t131 ^ 2;
	t161 = t121 * t126;
	t122 = sin(t124);
	t120 = t122 ^ 2;
	t132 = cos(pkin(9));
	t128 = 0.1e1 / t132 ^ 2;
	t116 = t120 * t126 * t128 + 0.1e1;
	t114 = 0.1e1 / t116;
	t160 = (t114 - 0.1e1) * t131;
	t133 = sin(qJ(5));
	t134 = cos(qJ(5));
	t147 = t132 * t134;
	t110 = t122 * t133 + t123 * t147;
	t151 = t122 * t131;
	t113 = atan2(-t151, -t132);
	t111 = sin(t113);
	t112 = cos(t113);
	t99 = -t111 * t151 - t112 * t132;
	t96 = 0.1e1 / t99;
	t104 = 0.1e1 / t110;
	t127 = 0.1e1 / t132;
	t105 = 0.1e1 / t110 ^ 2;
	t97 = 0.1e1 / t99 ^ 2;
	t148 = t132 * t133;
	t109 = -t122 * t134 + t123 * t148;
	t103 = t109 ^ 2;
	t102 = t103 * t105 + 0.1e1;
	t130 = qJD(1) + qJD(2);
	t143 = t122 * t148 + t123 * t134;
	t94 = -t110 * qJD(5) + t143 * t130;
	t154 = t94 * t105;
	t108 = -t122 * t147 + t123 * t133;
	t95 = -t109 * qJD(5) + t108 * t130;
	t157 = t104 * t105 * t95;
	t158 = 0.1e1 / t102 ^ 2 * (-t103 * t157 - t109 * t154);
	t156 = t122 * t97;
	t155 = t123 * t97;
	t153 = t108 * t109;
	t152 = t122 * t130;
	t149 = t126 * t127;
	t146 = t114 * t149;
	t115 = 0.1e1 / t116 ^ 2;
	t145 = t115 * t131 * t161;
	t89 = (-t112 * t122 * t146 + t111 * t160) * t123;
	t129 = t127 * t128;
	t100 = 0.1e1 / t102;
	t98 = t96 * t97;
	t93 = t97 * t161 + 0.1e1;
	t90 = (-t114 * t127 * t131 - 0.2e1 * t129 * t145) * t152;
	t88 = t130 * t89;
	t85 = 0.2e1 * (t104 * t143 + t105 * t153) * t158 + ((t108 * qJD(5) - t109 * t130) * t104 + 0.2e1 * t153 * t157 + (t143 * t95 - (t143 * qJD(5) - t110 * t130) * t109 + t108 * t94) * t105) * t100;
	t84 = (0.2e1 * (t122 * t96 + t89 * t155) / t93 ^ 2 * (-t121 * t88 * t98 - t152 * t155) * t126 + ((0.2e1 * t123 * t89 * t98 + t156) * t88 + (t89 * t156 + (-t96 + (t128 * t145 + t160) * t111 * t156 - (t120 * t146 + (-0.2e1 * t146 + (0.2e1 * t120 * t126 ^ 2 * t129 + t149) * t115) * t121) * t97 * t112) * t123) * t130) / t93) * t131;
	t1 = [t90, t90, 0, 0, 0; t84, t84, 0, 0, 0; t85, t85, 0, 0, -0.2e1 * t158 + 0.2e1 * (-t100 * t154 + (-t100 * t157 - t105 * t158) * t109) * t109;];
	JaD_rot = t1;
end