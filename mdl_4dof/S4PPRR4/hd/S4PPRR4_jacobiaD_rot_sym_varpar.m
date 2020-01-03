% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4PPRR4
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
%   Wie in S4PPRR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4PPRR4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PPRR4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_jacobiaD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:18:49
	% EndTime: 2019-12-31 16:18:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:18:49
	% EndTime: 2019-12-31 16:18:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:18:49
	% EndTime: 2019-12-31 16:18:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:18:49
	% EndTime: 2019-12-31 16:18:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:18:49
	% EndTime: 2019-12-31 16:18:49
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (973->43), mult. (1166->102), div. (229->12), fcn. (1315->9), ass. (0->60)
	t104 = sin(pkin(6));
	t103 = pkin(7) + qJ(3);
	t100 = cos(t103);
	t99 = sin(t103);
	t95 = t99 ^ 2;
	t130 = t95 / t100 ^ 2;
	t120 = 0.1e1 + t130;
	t143 = t104 * t120;
	t142 = qJD(3) * t99;
	t106 = sin(qJ(4));
	t107 = cos(qJ(4));
	t105 = cos(pkin(6));
	t125 = t100 * t105;
	t84 = t104 * t106 + t107 * t125;
	t128 = t104 * t99;
	t87 = atan2(-t128, -t100);
	t85 = sin(t87);
	t86 = cos(t87);
	t70 = -t86 * t100 - t85 * t128;
	t67 = 0.1e1 / t70;
	t80 = 0.1e1 / t84;
	t101 = t104 ^ 2;
	t90 = t101 * t130 + 0.1e1;
	t88 = 0.1e1 / t90;
	t72 = t88 * t143;
	t117 = t104 * t72 - 0.1e1;
	t127 = t104 - t72;
	t129 = t100 * t85;
	t131 = t86 * t99;
	t61 = -t117 * t131 - t127 * t129;
	t140 = 0.2e1 * t61;
	t68 = 0.1e1 / t70 ^ 2;
	t81 = 0.1e1 / t84 ^ 2;
	t102 = t105 ^ 2;
	t123 = qJD(3) * t100;
	t135 = t68 * t99;
	t121 = t86 * t128;
	t71 = qJD(3) * t72;
	t60 = (-t121 + t129) * t71 + (-t104 * t129 + t131) * qJD(3);
	t138 = t60 * t67 * t68;
	t66 = t102 * t95 * t68 + 0.1e1;
	t139 = (t123 * t135 - t95 * t138) * t102 / t66 ^ 2;
	t83 = -t104 * t107 + t106 * t125;
	t132 = t81 * t83;
	t119 = t105 * t142;
	t126 = qJD(4) * t83;
	t77 = -t107 * t119 - t126;
	t133 = t77 * t80 * t81;
	t79 = t83 ^ 2;
	t75 = t79 * t81 + 0.1e1;
	t76 = -t84 * qJD(4) + t106 * t119;
	t137 = (-t76 * t132 - t79 * t133) / t75 ^ 2;
	t64 = 0.1e1 / t66;
	t136 = t64 * t68;
	t134 = t76 * t81;
	t122 = -0.2e1 * t137;
	t116 = -t106 * t80 + t107 * t132;
	t73 = 0.1e1 / t75;
	t62 = 0.2e1 * (t88 - t120 / t90 ^ 2 * t101) / t100 * t142 * t143;
	t1 = [0, 0, t62, 0; 0, 0, ((-0.2e1 * t67 * t139 + (-qJD(3) * t61 - t60) * t136) * t100 + (t68 * t139 * t140 + (t62 * t68 * t121 + t138 * t140 - qJD(3) * t67 - (t127 * qJD(3) + t117 * t71) * t85 * t135) * t64 - (t62 * t85 + (qJD(3) + (-0.2e1 * t104 + t72) * t71) * t86) * t100 * t136) * t99) * t105, 0; 0, 0, (t116 * t99 * t122 + (t116 * t123 + ((t77 - t126) * t81 * t106 + (-qJD(4) * t80 - 0.2e1 * t83 * t133 - t134) * t107) * t99) * t73) * t105, t122 + 0.2e1 * (-t73 * t134 + (-t73 * t133 - t81 * t137) * t83) * t83;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,4);
end