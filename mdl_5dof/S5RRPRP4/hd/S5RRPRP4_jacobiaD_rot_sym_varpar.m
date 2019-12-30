% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPRP4
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
%   Wie in S5RRPRP4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRPRP4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_jacobiaD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:42:25
	% EndTime: 2019-12-29 18:42:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:42:25
	% EndTime: 2019-12-29 18:42:25
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:42:25
	% EndTime: 2019-12-29 18:42:25
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (74->0), mult. (74->0), div. (30->0), fcn. (44->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:42:20
	% EndTime: 2019-12-29 18:42:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:42:25
	% EndTime: 2019-12-29 18:42:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:42:08
	% EndTime: 2019-12-29 18:42:09
	% DurationCPUTime: 1.38s
	% Computational Cost: add. (1950->70), mult. (2645->163), div. (674->13), fcn. (3179->7), ass. (0->77)
	t112 = qJ(1) + qJ(2);
	t105 = cos(t112);
	t104 = sin(t112);
	t99 = 0.1e1 / t104;
	t153 = t105 * t99;
	t113 = sin(qJ(4));
	t108 = 0.1e1 / t113 ^ 2;
	t114 = cos(qJ(4));
	t111 = t114 ^ 2;
	t142 = t108 * t111;
	t128 = 0.1e1 + t142;
	t103 = t105 ^ 2;
	t97 = t103 * t142 + 0.1e1;
	t95 = 0.1e1 / t97;
	t125 = t128 * t95;
	t107 = 0.1e1 / t113;
	t131 = t105 * t142;
	t139 = qJD(4) * t113;
	t106 = qJD(1) + qJD(2);
	t145 = t106 * t114;
	t78 = ((t104 * t145 + t105 * t139) * t107 + qJD(4) * t131) * t95;
	t88 = t105 * t125;
	t164 = qJD(4) * t88 + t78;
	t146 = t105 * t114;
	t94 = atan2(-t146, t113);
	t91 = sin(t94);
	t92 = cos(t94);
	t87 = t92 * t113 - t146 * t91;
	t84 = 0.1e1 / t87;
	t85 = 0.1e1 / t87 ^ 2;
	t162 = t95 - 0.1e1;
	t132 = t85 * t139;
	t147 = t105 * t106;
	t135 = t85 * t147;
	t140 = qJD(4) * t105;
	t148 = t104 * t106;
	t152 = t113 * t91;
	t154 = t105 * t78;
	t73 = (-t78 + t140) * t152 + (t91 * t148 + (qJD(4) - t154) * t92) * t114;
	t160 = t73 * t84 * t85;
	t98 = t104 ^ 2;
	t81 = t111 * t85 * t98 + 0.1e1;
	t161 = (-t98 * t114 * t132 + (t104 * t135 - t160 * t98) * t111) / t81 ^ 2;
	t79 = 0.1e1 / t81;
	t159 = t79 * t85;
	t109 = t107 * t108;
	t123 = t106 * (-0.1e1 / t98 * t103 - 0.1e1) * t153;
	t149 = 0.1e1 / t104 ^ 2 * t103;
	t93 = t108 * t149 + 0.1e1;
	t158 = (-qJD(4) * t109 * t114 * t149 + t108 * t123) / t93 ^ 2;
	t110 = t114 * t111;
	t143 = t107 * t114;
	t122 = qJD(4) * (-t109 * t110 - t143);
	t157 = (t103 * t122 - t131 * t148) / t97 ^ 2;
	t156 = t84 * t79;
	t155 = t104 * t85;
	t151 = t92 * t114;
	t144 = t107 * t111;
	t141 = t108 * t114;
	t138 = -0.2e1 * t160;
	t137 = t84 * t161;
	t136 = t114 * t157;
	t134 = t95 * t144;
	t133 = t114 * t162;
	t130 = -0.2e1 * t85 * t161;
	t129 = -0.1e1 - t149;
	t127 = t91 * t133;
	t126 = t105 * t134;
	t124 = t129 * t141;
	t89 = 0.1e1 / t93;
	t77 = (-t126 * t92 - t127) * t104;
	t76 = t95 * t143 * t147 + (-qJD(4) * t125 - 0.2e1 * t107 * t136) * t104;
	t75 = -t88 * t152 + t151 + (-t151 * t88 + t152) * t105;
	t74 = -t125 * t148 + 0.2e1 * (t122 * t95 - t128 * t157) * t105;
	t72 = t89 * qJD(4) * t124 + 0.2e1 * (t123 * t89 + t129 * t158) * t107;
	t70 = (t139 * t156 + (0.2e1 * t137 + (t106 * t77 + t73) * t159) * t114) * t105 + (t77 * t130 * t114 + (-t77 * t132 + (t77 * t138 + ((t126 * t78 + t139 * t162 + 0.2e1 * t136) * t91 + (-t78 * t133 + (0.2e1 * t144 * t157 + (t108 * t110 + 0.2e1 * t114) * t95 * qJD(4)) * t105) * t92) * t155) * t114 + (t84 + ((-t103 + t98) * t92 * t134 - t105 * t127) * t85) * t145) * t79) * t104;
	t1 = [t76, t76, 0, t74, 0; t70, t70, 0, (t147 * t156 + (-0.2e1 * t137 + (-qJD(4) * t75 - t73) * t159) * t104) * t113 + (t75 * t104 * t130 + (t104 * qJD(4) * t84 + (t104 * t138 + t135) * t75 + (((-t74 - t148) * t91 + (t164 * t105 - t78 * t88 - qJD(4)) * t92) * t113 + ((-t105 * t74 + t148 * t88) * t92 + (t154 * t88 + t140 - t164) * t91) * t114) * t155) * t79) * t114, 0; t72, t72, 0, -0.2e1 * t141 * t153 * t158 + (t106 * t124 + (-0.2e1 * t109 * t111 - t107) * t99 * t140) * t89, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end