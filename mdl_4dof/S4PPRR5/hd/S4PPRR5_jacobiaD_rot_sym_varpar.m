% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4PPRR5
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
%   Wie in S4PPRR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 11:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4PPRR5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PPRR5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_jacobiaD_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 11:59:16
	% EndTime: 2019-12-29 11:59:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 11:59:21
	% EndTime: 2019-12-29 11:59:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 11:59:16
	% EndTime: 2019-12-29 11:59:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 11:59:21
	% EndTime: 2019-12-29 11:59:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 11:59:21
	% EndTime: 2019-12-29 11:59:22
	% DurationCPUTime: 0.61s
	% Computational Cost: add. (355->40), mult. (1166->102), div. (229->12), fcn. (1315->9), ass. (0->55)
	t97 = cos(qJ(3));
	t91 = t97 ^ 2;
	t95 = sin(qJ(3));
	t119 = 0.1e1 / t95 ^ 2 * t91;
	t109 = 0.1e1 + t119;
	t93 = cos(pkin(6));
	t127 = t109 * t93;
	t115 = qJD(3) * t97;
	t117 = t93 * t97;
	t80 = atan2(-t117, t95);
	t78 = sin(t80);
	t79 = cos(t80);
	t64 = -t78 * t117 + t79 * t95;
	t61 = 0.1e1 / t64;
	t92 = sin(pkin(6));
	t118 = t92 * t95;
	t94 = sin(qJ(4));
	t96 = cos(qJ(4));
	t77 = t96 * t118 + t93 * t94;
	t73 = 0.1e1 / t77;
	t62 = 0.1e1 / t64 ^ 2;
	t74 = 0.1e1 / t77 ^ 2;
	t120 = t78 * t95;
	t106 = t93 * t120 + t79 * t97;
	t112 = t79 * t117;
	t107 = -t112 - t120;
	t86 = t93 ^ 2;
	t83 = t86 * t119 + 0.1e1;
	t81 = 0.1e1 / t83;
	t67 = t81 * t127;
	t60 = qJD(3) * t67;
	t53 = t106 * qJD(3) + t107 * t60;
	t125 = t53 * t61 * t62;
	t76 = t94 * t118 - t93 * t96;
	t121 = t74 * t76;
	t111 = t92 * t115;
	t114 = qJD(4) * t76;
	t70 = t96 * t111 - t114;
	t122 = t70 * t73 * t74;
	t72 = t76 ^ 2;
	t68 = t72 * t74 + 0.1e1;
	t69 = t77 * qJD(4) + t94 * t111;
	t124 = (t69 * t121 - t72 * t122) / t68 ^ 2;
	t123 = t62 * t95;
	t116 = -t67 + t93;
	t85 = t92 ^ 2;
	t59 = t85 * t91 * t62 + 0.1e1;
	t113 = -0.2e1 * (-t115 * t123 - t91 * t125) * t85 / t59 ^ 2;
	t110 = t67 * t93 - 0.1e1;
	t108 = t96 * t121 - t73 * t94;
	t65 = 0.1e1 / t68;
	t57 = 0.1e1 / t59;
	t56 = -0.2e1 * (t81 - t109 / t83 ^ 2 * t86) / t95 * t115 * t127;
	t54 = t107 * t67 + t106;
	t1 = [0, 0, t56, 0; 0, 0, ((t61 * t113 + (-qJD(3) * t54 - t53) * t62 * t57) * t95 + (t54 * t62 * t113 + (-0.2e1 * t54 * t125 + qJD(3) * t61 + (-t56 * t78 + (t110 * qJD(3) + t116 * t60) * t79) * t123 + (-t56 * t112 + (t116 * qJD(3) + t110 * t60) * t97 * t78) * t62) * t57) * t97) * t92, 0; 0, 0, (0.2e1 * t108 * t97 * t124 + (t108 * t95 * qJD(3) + ((qJD(4) * t73 + 0.2e1 * t76 * t122) * t96 + (-t69 * t96 + (-t70 + t114) * t94) * t74) * t97) * t65) * t92, -0.2e1 * t124 + 0.2e1 * (t65 * t69 * t74 + (-t65 * t122 - t74 * t124) * t76) * t76;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,4);
end