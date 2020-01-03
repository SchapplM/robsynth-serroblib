% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRPP1
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
%   Wie in S5RPRPP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPRPP1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPP1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:09:52
	% EndTime: 2019-12-31 18:09:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:09:52
	% EndTime: 2019-12-31 18:09:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:09:52
	% EndTime: 2019-12-31 18:09:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:09:52
	% EndTime: 2019-12-31 18:09:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:09:52
	% EndTime: 2019-12-31 18:09:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:09:52
	% EndTime: 2019-12-31 18:09:53
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (2575->73), mult. (1835->159), div. (470->13), fcn. (2177->7), ass. (0->73)
	t99 = qJ(1) + pkin(7);
	t95 = sin(t99);
	t128 = qJD(1) * t95;
	t148 = 0.2e1 * t95;
	t86 = t95 ^ 2;
	t97 = cos(t99);
	t90 = t97 ^ 2;
	t91 = 0.1e1 / t97;
	t146 = (t86 / t90 + 0.1e1) * t91 * t128;
	t98 = qJ(3) + pkin(8);
	t94 = sin(t98);
	t129 = t95 * t94;
	t96 = cos(t98);
	t76 = atan2(-t129, -t96);
	t74 = sin(t76);
	t120 = t74 * t129;
	t75 = cos(t76);
	t70 = -t75 * t96 - t120;
	t67 = 0.1e1 / t70;
	t87 = 0.1e1 / t96;
	t145 = -0.2e1 * t94;
	t68 = 0.1e1 / t70 ^ 2;
	t88 = 0.1e1 / t96 ^ 2;
	t84 = t94 ^ 2;
	t133 = t84 * t88;
	t81 = t133 * t86 + 0.1e1;
	t77 = 0.1e1 / t81;
	t144 = t77 - 0.1e1;
	t127 = qJD(1) * t97;
	t118 = t94 * t127;
	t126 = qJD(3) * t95;
	t135 = t75 * t94;
	t125 = qJD(3) * t96;
	t63 = (-(-t125 * t95 - t118) * t87 + t126 * t133) * t77;
	t59 = (-t63 * t95 + qJD(3)) * t135 + (-t118 + (t63 - t126) * t96) * t74;
	t143 = t59 * t67 * t68;
	t142 = t63 * t74;
	t141 = t63 * t94;
	t140 = t68 * t94;
	t139 = t68 * t97;
	t130 = t87 * t94;
	t83 = t94 * t84;
	t89 = t87 * t88;
	t107 = qJD(3) * (t83 * t89 + t130);
	t111 = t84 * t95 * t127;
	t138 = (t107 * t86 + t111 * t88) / t81 ^ 2;
	t117 = 0.1e1 + t133;
	t73 = t117 * t95 * t77;
	t137 = t73 * t95;
	t136 = t74 * t96;
	t134 = t84 * t87;
	t132 = t84 * t90;
	t131 = t86 / t97 ^ 2;
	t66 = t132 * t68 + 0.1e1;
	t124 = 0.2e1 * (-t132 * t143 + (t125 * t90 * t94 - t111) * t68) / t66 ^ 2;
	t123 = 0.2e1 * t143;
	t82 = t131 * t88 + 0.1e1;
	t122 = 0.2e1 * (qJD(3) * t131 * t89 * t94 + t146 * t88) / t82 ^ 2;
	t121 = t94 * t139;
	t119 = t77 * t134;
	t116 = 0.1e1 + t131;
	t115 = t94 * t124;
	t114 = t138 * t145;
	t113 = t138 * t148;
	t112 = t95 * t119;
	t110 = t117 * t97;
	t108 = t116 * t94 * t88;
	t79 = 0.1e1 / t82;
	t64 = 0.1e1 / t66;
	t62 = (t144 * t74 * t94 - t112 * t75) * t97;
	t61 = -t95 * t136 + t135 + (-t129 * t75 + t136) * t73;
	t60 = -t117 * t113 + (qJD(1) * t110 + t107 * t148) * t77;
	t1 = [t97 * t87 * t114 + (qJD(3) * t110 - t128 * t130) * t77, 0, t60, 0, 0; (t67 * t115 + (-t67 * t125 + (qJD(1) * t62 + t59) * t140) * t64) * t95 + (t68 * t115 * t62 + (-((t112 * t63 + t125 * t144 + t114) * t74 + (t113 * t134 - t141 + (t141 + (-t83 * t88 + t145) * t126) * t77) * t75) * t121 + (t123 * t94 - t125 * t68) * t62 + (-t67 + ((-t86 + t90) * t75 * t119 + t144 * t120) * t68) * t94 * qJD(1)) * t64) * t97, 0, (t140 * t61 - t67 * t96) * t97 * t124 + ((-t67 * t128 + (-qJD(3) * t61 - t59) * t139) * t96 + (-t97 * qJD(3) * t67 - (-t60 * t75 * t95 + t74 * t126 + t137 * t142 - t142 + (-qJD(3) * t74 - t127 * t75) * t73) * t121 + (t123 * t97 + t128 * t68) * t61 - ((t60 - t127) * t74 + ((0.1e1 - t137) * qJD(3) + (t73 - t95) * t63) * t75) * t96 * t139) * t94) * t64, 0, 0; t116 * t87 * t122 + (-qJD(3) * t108 - 0.2e1 * t146 * t87) * t79, 0, t91 * t88 * t122 * t129 + ((-0.2e1 * t84 * t89 - t87) * t91 * t126 - qJD(1) * t108) * t79, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end