% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRPP1
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
%   Wie in S5RRRPP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRPP1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:33:47
	% EndTime: 2019-12-29 19:33:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:33:37
	% EndTime: 2019-12-29 19:33:37
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:33:42
	% EndTime: 2019-12-29 19:33:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (74->0), mult. (74->0), div. (30->0), fcn. (44->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:33:32
	% EndTime: 2019-12-29 19:33:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:33:42
	% EndTime: 2019-12-29 19:33:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:33:42
	% EndTime: 2019-12-29 19:33:43
	% DurationCPUTime: 1.51s
	% Computational Cost: add. (3852->72), mult. (2645->160), div. (674->13), fcn. (3179->7), ass. (0->74)
	t117 = qJ(1) + qJ(2);
	t113 = sin(t117);
	t115 = qJD(1) + qJD(2);
	t145 = t113 * t115;
	t116 = qJ(3) + pkin(8);
	t111 = sin(t116);
	t101 = t111 ^ 2;
	t112 = cos(t116);
	t103 = 0.1e1 / t112 ^ 2;
	t148 = t101 * t103;
	t132 = 0.1e1 + t148;
	t106 = t113 ^ 2;
	t98 = t106 * t148 + 0.1e1;
	t94 = 0.1e1 / t98;
	t128 = t132 * t94;
	t90 = t113 * t128;
	t165 = t113 * t90 - 0.1e1;
	t114 = cos(t117);
	t107 = t114 ^ 2;
	t108 = 0.1e1 / t114;
	t163 = (t106 / t107 + 0.1e1) * t108 * t145;
	t146 = t113 * t111;
	t93 = atan2(-t146, -t112);
	t91 = sin(t93);
	t137 = t91 * t146;
	t92 = cos(t93);
	t89 = -t112 * t92 - t137;
	t86 = 0.1e1 / t89;
	t102 = 0.1e1 / t112;
	t87 = 0.1e1 / t89 ^ 2;
	t162 = 0.2e1 * t111;
	t161 = t94 - 0.1e1;
	t144 = t114 * t115;
	t134 = t113 * t144;
	t143 = qJD(3) * t112;
	t135 = t87 * t143;
	t142 = qJD(3) * t113;
	t152 = t112 * t91;
	t80 = (-(-t111 * t144 - t112 * t142) * t102 + t142 * t148) * t94;
	t74 = (t80 - t142) * t152 + (-t91 * t144 + (-t113 * t80 + qJD(3)) * t92) * t111;
	t159 = t74 * t86 * t87;
	t154 = t101 * t87;
	t83 = t107 * t154 + 0.1e1;
	t160 = (t107 * t111 * t135 + (-t107 * t159 - t134 * t87) * t101) / t83 ^ 2;
	t81 = 0.1e1 / t83;
	t157 = t81 * t87;
	t100 = t111 * t101;
	t104 = t102 * t103;
	t125 = qJD(3) * (t100 * t104 + t102 * t111);
	t156 = (t106 * t125 + t134 * t148) / t98 ^ 2;
	t155 = t86 * t81;
	t153 = t102 * t94;
	t150 = t114 * t87;
	t149 = qJD(3) * t90;
	t147 = t106 / t114 ^ 2;
	t141 = 0.2e1 * t159;
	t99 = t103 * t147 + 0.1e1;
	t140 = 0.2e1 * (qJD(3) * t104 * t111 * t147 + t103 * t163) / t99 ^ 2;
	t139 = t86 * t160;
	t138 = t113 * t153;
	t136 = t161 * t111;
	t133 = 0.2e1 * t87 * t160;
	t131 = 0.1e1 + t147;
	t130 = -0.2e1 * t102 * t156;
	t129 = t101 * t138;
	t126 = t103 * t111 * t131;
	t96 = 0.1e1 / t99;
	t79 = (-t129 * t92 + t136 * t91) * t114;
	t78 = -t111 * t115 * t138 + (qJD(3) * t128 + t111 * t130) * t114;
	t77 = (-t113 + t90) * t152 - t165 * t92 * t111;
	t76 = t128 * t144 + 0.2e1 * (t125 * t94 - t132 * t156) * t113;
	t75 = -t96 * qJD(3) * t126 + (t131 * t140 - 0.2e1 * t96 * t163) * t102;
	t72 = (-t143 * t155 + (0.2e1 * t139 + (t115 * t79 + t74) * t157) * t111) * t113 + (-t79 * t81 * t135 + (t79 * t133 + (t79 * t141 + ((-t129 * t80 - t143 * t161 + t156 * t162) * t91 + (-t80 * t136 + (t101 * t130 + (t100 * t103 + t162) * t94 * qJD(3)) * t113) * t92) * t150 + (-t86 + t161 * t87 * t137 - (t106 - t107) * t92 * t153 * t154) * t115) * t81) * t111) * t114;
	t1 = [t78, t78, t76, 0, 0; t72, t72, (-t145 * t155 + (-0.2e1 * t139 + (-qJD(3) * t77 - t74) * t157) * t114) * t112 + (t77 * t114 * t133 + (-t114 * qJD(3) * t86 + (t114 * t141 + t145 * t87) * t77 + (-((-t113 * t76 - t144 * t90) * t92 + (t165 * t80 + t142 - t149) * t91) * t111 - ((t76 - t144) * t91 + (t80 * t90 + qJD(3) + (-t80 - t149) * t113) * t92) * t112) * t150) * t81) * t111, 0, 0; t75, t75, t103 * t108 * t140 * t146 + (-t115 * t126 + (-0.2e1 * t101 * t104 - t102) * t108 * t142) * t96, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end