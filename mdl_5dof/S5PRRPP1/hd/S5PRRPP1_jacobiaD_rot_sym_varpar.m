% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRPP1
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
%   Wie in S5PRRPP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PRRPP1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPP1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:07:41
	% EndTime: 2019-12-05 16:07:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:07:41
	% EndTime: 2019-12-05 16:07:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:07:41
	% EndTime: 2019-12-05 16:07:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:07:41
	% EndTime: 2019-12-05 16:07:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:07:41
	% EndTime: 2019-12-05 16:07:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:07:41
	% EndTime: 2019-12-05 16:07:42
	% DurationCPUTime: 0.60s
	% Computational Cost: add. (2575->73), mult. (1835->160), div. (470->13), fcn. (2177->7), ass. (0->73)
	t98 = pkin(7) + qJ(2);
	t94 = sin(t98);
	t129 = qJD(2) * t94;
	t148 = 0.2e1 * t94;
	t84 = t94 ^ 2;
	t96 = cos(t98);
	t87 = t96 ^ 2;
	t88 = 0.1e1 / t96;
	t146 = (t84 / t87 + 0.1e1) * t88 * t129;
	t99 = qJ(3) + pkin(8);
	t95 = sin(t99);
	t130 = t94 * t95;
	t97 = cos(t99);
	t76 = atan2(-t130, -t97);
	t74 = sin(t76);
	t120 = t74 * t130;
	t75 = cos(t76);
	t70 = -t75 * t97 - t120;
	t67 = 0.1e1 / t70;
	t91 = 0.1e1 / t97;
	t145 = -0.2e1 * t95;
	t68 = 0.1e1 / t70 ^ 2;
	t92 = 0.1e1 / t97 ^ 2;
	t86 = t95 ^ 2;
	t131 = t86 * t92;
	t81 = t131 * t84 + 0.1e1;
	t77 = 0.1e1 / t81;
	t144 = t77 - 0.1e1;
	t127 = qJD(2) * t96;
	t118 = t95 * t127;
	t126 = qJD(3) * t94;
	t135 = t75 * t95;
	t125 = qJD(3) * t97;
	t63 = (-(-t125 * t94 - t118) * t91 + t126 * t131) * t77;
	t59 = (-t63 * t94 + qJD(3)) * t135 + (-t118 + (t63 - t126) * t97) * t74;
	t143 = t59 * t67 * t68;
	t142 = t63 * t74;
	t141 = t63 * t95;
	t140 = t68 * t95;
	t139 = t68 * t96;
	t85 = t95 * t86;
	t93 = t91 * t92;
	t107 = qJD(3) * (t85 * t93 + t91 * t95);
	t111 = t86 * t94 * t127;
	t138 = (t107 * t84 + t111 * t92) / t81 ^ 2;
	t116 = 0.1e1 + t131;
	t73 = t116 * t94 * t77;
	t137 = t73 * t94;
	t136 = t74 * t97;
	t134 = t84 / t96 ^ 2;
	t133 = t86 * t87;
	t132 = t86 * t91;
	t128 = qJD(2) * t95;
	t66 = t133 * t68 + 0.1e1;
	t124 = 0.2e1 * (-t133 * t143 + (t125 * t87 * t95 - t111) * t68) / t66 ^ 2;
	t123 = 0.2e1 * t143;
	t82 = t134 * t92 + 0.1e1;
	t122 = 0.2e1 * (qJD(3) * t134 * t93 * t95 + t146 * t92) / t82 ^ 2;
	t121 = t95 * t139;
	t119 = t77 * t132;
	t117 = 0.1e1 + t134;
	t115 = t95 * t124;
	t114 = t138 * t148;
	t113 = t138 * t145;
	t112 = t94 * t119;
	t110 = t116 * t96;
	t108 = t117 * t95 * t92;
	t79 = 0.1e1 / t82;
	t64 = 0.1e1 / t66;
	t62 = (t144 * t74 * t95 - t112 * t75) * t96;
	t61 = -t94 * t136 + t135 + (-t130 * t75 + t136) * t73;
	t60 = -t116 * t114 + (qJD(2) * t110 + t107 * t148) * t77;
	t1 = [0, t96 * t91 * t113 + (-t128 * t91 * t94 + qJD(3) * t110) * t77, t60, 0, 0; 0, (t67 * t115 + (-t67 * t125 + (qJD(2) * t62 + t59) * t140) * t64) * t94 + (t68 * t115 * t62 + (-((t112 * t63 + t125 * t144 + t113) * t74 + (t114 * t132 - t141 + (t141 + (-t85 * t92 + t145) * t126) * t77) * t75) * t121 + (t123 * t95 - t125 * t68) * t62 + (-t67 + ((-t84 + t87) * t75 * t119 + t144 * t120) * t68) * t128) * t64) * t96, (t140 * t61 - t67 * t97) * t96 * t124 + ((-t67 * t129 + (-qJD(3) * t61 - t59) * t139) * t97 + (-t96 * qJD(3) * t67 - (-t60 * t75 * t94 + t74 * t126 + t137 * t142 - t142 + (-qJD(3) * t74 - t127 * t75) * t73) * t121 + (t123 * t96 + t129 * t68) * t61 - ((t60 - t127) * t74 + ((0.1e1 - t137) * qJD(3) + (t73 - t94) * t63) * t75) * t97 * t139) * t95) * t64, 0, 0; 0, t117 * t91 * t122 + (-qJD(3) * t108 - 0.2e1 * t146 * t91) * t79, t88 * t92 * t122 * t130 + ((-0.2e1 * t86 * t93 - t91) * t88 * t126 - qJD(2) * t108) * t79, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end