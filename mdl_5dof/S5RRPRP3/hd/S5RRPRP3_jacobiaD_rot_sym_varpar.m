% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPRP3
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
%   Wie in S5RRPRP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRPRP3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:39:33
	% EndTime: 2019-12-29 18:39:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:39:34
	% EndTime: 2019-12-29 18:39:34
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:39:39
	% EndTime: 2019-12-29 18:39:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (74->0), mult. (74->0), div. (30->0), fcn. (44->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:39:39
	% EndTime: 2019-12-29 18:39:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:39:34
	% EndTime: 2019-12-29 18:39:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:39:39
	% EndTime: 2019-12-29 18:39:40
	% DurationCPUTime: 1.40s
	% Computational Cost: add. (3852->72), mult. (2645->159), div. (674->13), fcn. (3179->7), ass. (0->74)
	t115 = qJ(1) + qJ(2);
	t111 = sin(t115);
	t114 = qJD(1) + qJD(2);
	t143 = t111 * t114;
	t113 = pkin(8) + qJ(4);
	t110 = cos(t113);
	t101 = 0.1e1 / t110 ^ 2;
	t109 = sin(t113);
	t99 = t109 ^ 2;
	t150 = t101 * t99;
	t132 = 0.1e1 + t150;
	t104 = t111 ^ 2;
	t96 = t104 * t150 + 0.1e1;
	t92 = 0.1e1 / t96;
	t126 = t132 * t92;
	t88 = t111 * t126;
	t163 = t111 * t88 - 0.1e1;
	t112 = cos(t115);
	t105 = t112 ^ 2;
	t106 = 0.1e1 / t112;
	t161 = (t104 / t105 + 0.1e1) * t106 * t143;
	t144 = t111 * t109;
	t91 = atan2(-t144, -t110);
	t89 = sin(t91);
	t135 = t89 * t144;
	t90 = cos(t91);
	t87 = -t110 * t90 - t135;
	t84 = 0.1e1 / t87;
	t100 = 0.1e1 / t110;
	t85 = 0.1e1 / t87 ^ 2;
	t160 = 0.2e1 * t109;
	t159 = t92 - 0.1e1;
	t142 = t112 * t114;
	t127 = t111 * t99 * t142;
	t141 = qJD(4) * t110;
	t133 = t85 * t141;
	t140 = qJD(4) * t111;
	t149 = t110 * t89;
	t78 = (-(-t109 * t142 - t110 * t140) * t100 + t140 * t150) * t92;
	t72 = (t78 - t140) * t149 + (-t89 * t142 + (-t111 * t78 + qJD(4)) * t90) * t109;
	t157 = t72 * t84 * t85;
	t152 = t85 * t99;
	t81 = t105 * t152 + 0.1e1;
	t158 = (-t85 * t127 + (t109 * t133 - t157 * t99) * t105) / t81 ^ 2;
	t79 = 0.1e1 / t81;
	t155 = t79 * t85;
	t102 = t100 * t101;
	t98 = t109 * t99;
	t123 = qJD(4) * (t100 * t109 + t102 * t98);
	t154 = (t101 * t127 + t104 * t123) / t96 ^ 2;
	t153 = t84 * t79;
	t151 = t100 * t92;
	t147 = t112 * t85;
	t146 = qJD(4) * t88;
	t145 = t104 / t112 ^ 2;
	t139 = 0.2e1 * t157;
	t97 = t101 * t145 + 0.1e1;
	t138 = 0.2e1 * (qJD(4) * t102 * t109 * t145 + t101 * t161) / t97 ^ 2;
	t137 = t84 * t158;
	t136 = t111 * t151;
	t134 = t159 * t109;
	t131 = 0.2e1 * t85 * t158;
	t130 = 0.1e1 + t145;
	t129 = -0.2e1 * t100 * t154;
	t128 = t99 * t136;
	t124 = t101 * t109 * t130;
	t94 = 0.1e1 / t97;
	t77 = (-t128 * t90 + t134 * t89) * t112;
	t76 = -t109 * t114 * t136 + (qJD(4) * t126 + t109 * t129) * t112;
	t75 = (-t111 + t88) * t149 - t163 * t90 * t109;
	t74 = t126 * t142 + 0.2e1 * (t123 * t92 - t132 * t154) * t111;
	t73 = -t94 * qJD(4) * t124 + (t130 * t138 - 0.2e1 * t94 * t161) * t100;
	t70 = (-t141 * t153 + (0.2e1 * t137 + (t114 * t77 + t72) * t155) * t109) * t111 + (-t77 * t79 * t133 + (t77 * t131 + (t77 * t139 + ((-t78 * t128 - t159 * t141 + t154 * t160) * t89 + (-t78 * t134 + (t99 * t129 + (t101 * t98 + t160) * t92 * qJD(4)) * t111) * t90) * t147 + (-t84 + t159 * t85 * t135 - (t104 - t105) * t90 * t151 * t152) * t114) * t79) * t109) * t112;
	t1 = [t76, t76, 0, t74, 0; t70, t70, 0, (-t143 * t153 + (-0.2e1 * t137 + (-qJD(4) * t75 - t72) * t155) * t112) * t110 + (t75 * t112 * t131 + (-t112 * qJD(4) * t84 + (t112 * t139 + t143 * t85) * t75 + (-((-t111 * t74 - t142 * t88) * t90 + (t163 * t78 + t140 - t146) * t89) * t109 - ((t74 - t142) * t89 + (t78 * t88 + qJD(4) + (-t78 - t146) * t111) * t90) * t110) * t147) * t79) * t109, 0; t73, t73, 0, t101 * t106 * t138 * t144 + (-t114 * t124 + (-0.2e1 * t102 * t99 - t100) * t106 * t140) * t94, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end