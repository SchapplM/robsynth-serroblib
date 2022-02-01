% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRPR3
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
%   Wie in S5RPRPR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPRPR3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
JaD_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:21:20
	% EndTime: 2022-01-23 09:21:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:21:20
	% EndTime: 2022-01-23 09:21:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:21:20
	% EndTime: 2022-01-23 09:21:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:21:20
	% EndTime: 2022-01-23 09:21:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (118->0), mult. (74->0), div. (30->0), fcn. (44->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:21:20
	% EndTime: 2022-01-23 09:21:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:21:20
	% EndTime: 2022-01-23 09:21:20
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (1990->41), mult. (1512->108), div. (234->12), fcn. (1884->9), ass. (0->57)
	t122 = qJ(1) + pkin(8) + qJ(3);
	t121 = cos(t122);
	t119 = t121 ^ 2;
	t129 = sin(pkin(9));
	t124 = t129 ^ 2;
	t159 = t119 * t124;
	t120 = sin(t122);
	t118 = t120 ^ 2;
	t130 = cos(pkin(9));
	t126 = 0.1e1 / t130 ^ 2;
	t114 = t118 * t124 * t126 + 0.1e1;
	t112 = 0.1e1 / t114;
	t158 = (t112 - 0.1e1) * t129;
	t131 = sin(qJ(5));
	t132 = cos(qJ(5));
	t145 = t130 * t132;
	t108 = t120 * t131 + t121 * t145;
	t149 = t120 * t129;
	t111 = atan2(-t149, -t130);
	t109 = sin(t111);
	t110 = cos(t111);
	t97 = -t109 * t149 - t110 * t130;
	t94 = 0.1e1 / t97;
	t102 = 0.1e1 / t108;
	t125 = 0.1e1 / t130;
	t103 = 0.1e1 / t108 ^ 2;
	t95 = 0.1e1 / t97 ^ 2;
	t146 = t130 * t131;
	t107 = -t120 * t132 + t121 * t146;
	t101 = t107 ^ 2;
	t100 = t101 * t103 + 0.1e1;
	t128 = qJD(1) + qJD(3);
	t141 = t120 * t146 + t121 * t132;
	t92 = -t108 * qJD(5) + t141 * t128;
	t152 = t92 * t103;
	t106 = -t120 * t145 + t121 * t131;
	t93 = -t107 * qJD(5) + t106 * t128;
	t155 = t102 * t103 * t93;
	t156 = (-t101 * t155 - t107 * t152) / t100 ^ 2;
	t154 = t120 * t95;
	t153 = t121 * t95;
	t151 = t106 * t107;
	t150 = t120 * t128;
	t147 = t124 * t125;
	t144 = t112 * t147;
	t113 = 0.1e1 / t114 ^ 2;
	t143 = t113 * t129 * t159;
	t87 = (-t110 * t120 * t144 + t109 * t158) * t121;
	t127 = t125 * t126;
	t98 = 0.1e1 / t100;
	t96 = t94 * t95;
	t91 = t95 * t159 + 0.1e1;
	t88 = (-t112 * t125 * t129 - 0.2e1 * t127 * t143) * t150;
	t86 = t128 * t87;
	t83 = 0.2e1 * (t102 * t141 + t103 * t151) * t156 + ((t106 * qJD(5) - t107 * t128) * t102 + 0.2e1 * t151 * t155 + (t141 * t93 - (t141 * qJD(5) - t108 * t128) * t107 + t106 * t92) * t103) * t98;
	t82 = (0.2e1 * (t120 * t94 + t87 * t153) / t91 ^ 2 * (-t119 * t86 * t96 - t150 * t153) * t124 + ((0.2e1 * t121 * t87 * t96 + t154) * t86 + (t87 * t154 + (-t94 + (t126 * t143 + t158) * t109 * t154 - (t118 * t144 + (-0.2e1 * t144 + (0.2e1 * t118 * t124 ^ 2 * t127 + t147) * t113) * t119) * t95 * t110) * t121) * t128) / t91) * t129;
	t1 = [t88, 0, t88, 0, 0; t82, 0, t82, 0, 0; t83, 0, t83, 0, -0.2e1 * t156 + 0.2e1 * (-t98 * t152 + (-t103 * t156 - t98 * t155) * t107) * t107;];
	JaD_rot = t1;
end