% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4PRPR5
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S4PRPR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4PRPR5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRPR5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_jacobiaD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:23:27
	% EndTime: 2019-12-31 16:23:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:23:27
	% EndTime: 2019-12-31 16:23:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:23:27
	% EndTime: 2019-12-31 16:23:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:23:27
	% EndTime: 2019-12-31 16:23:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:23:28
	% EndTime: 2019-12-31 16:23:28
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (973->43), mult. (1166->102), div. (229->12), fcn. (1315->9), ass. (0->60)
	t105 = sin(pkin(6));
	t104 = qJ(2) + pkin(7);
	t101 = cos(t104);
	t100 = sin(t104);
	t96 = t100 ^ 2;
	t133 = t96 / t101 ^ 2;
	t121 = 0.1e1 + t133;
	t144 = t105 * t121;
	t143 = qJD(2) * t100;
	t107 = sin(qJ(4));
	t108 = cos(qJ(4));
	t106 = cos(pkin(6));
	t127 = t101 * t106;
	t85 = t105 * t107 + t108 * t127;
	t126 = t105 * t100;
	t88 = atan2(-t126, -t101);
	t86 = sin(t88);
	t87 = cos(t88);
	t71 = -t87 * t101 - t86 * t126;
	t68 = 0.1e1 / t71;
	t81 = 0.1e1 / t85;
	t102 = t105 ^ 2;
	t91 = t102 * t133 + 0.1e1;
	t89 = 0.1e1 / t91;
	t73 = t89 * t144;
	t119 = t105 * t73 - 0.1e1;
	t129 = t105 - t73;
	t130 = t87 * t100;
	t131 = t101 * t86;
	t62 = -t119 * t130 - t129 * t131;
	t141 = 0.2e1 * t62;
	t69 = 0.1e1 / t71 ^ 2;
	t82 = 0.1e1 / t85 ^ 2;
	t103 = t106 ^ 2;
	t124 = qJD(2) * t101;
	t132 = t100 * t69;
	t122 = t87 * t126;
	t72 = qJD(2) * t73;
	t61 = (-t122 + t131) * t72 + (-t105 * t131 + t130) * qJD(2);
	t139 = t61 * t68 * t69;
	t67 = t103 * t96 * t69 + 0.1e1;
	t140 = (t124 * t132 - t96 * t139) * t103 / t67 ^ 2;
	t84 = -t105 * t108 + t107 * t127;
	t134 = t82 * t84;
	t118 = t106 * t143;
	t128 = qJD(4) * t84;
	t78 = -t108 * t118 - t128;
	t135 = t78 * t81 * t82;
	t80 = t84 ^ 2;
	t76 = t80 * t82 + 0.1e1;
	t77 = -t85 * qJD(4) + t107 * t118;
	t138 = (-t77 * t134 - t80 * t135) / t76 ^ 2;
	t65 = 0.1e1 / t67;
	t137 = t65 * t69;
	t136 = t77 * t82;
	t123 = -0.2e1 * t138;
	t117 = -t107 * t81 + t108 * t134;
	t74 = 0.1e1 / t76;
	t63 = 0.2e1 * (t89 - t121 / t91 ^ 2 * t102) / t101 * t143 * t144;
	t1 = [0, t63, 0, 0; 0, ((-0.2e1 * t68 * t140 + (-qJD(2) * t62 - t61) * t137) * t101 + (t69 * t140 * t141 + (t63 * t69 * t122 + t139 * t141 - qJD(2) * t68 - (t129 * qJD(2) + t119 * t72) * t86 * t132) * t65 - (t63 * t86 + (qJD(2) + (-0.2e1 * t105 + t73) * t72) * t87) * t101 * t137) * t100) * t106, 0, 0; 0, (t117 * t74 * t124 + (t117 * t123 + ((t78 - t128) * t82 * t107 + (-qJD(4) * t81 - 0.2e1 * t84 * t135 - t136) * t108) * t74) * t100) * t106, 0, t123 + 0.2e1 * (-t74 * t136 + (-t74 * t135 - t82 * t138) * t84) * t84;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,4);
end