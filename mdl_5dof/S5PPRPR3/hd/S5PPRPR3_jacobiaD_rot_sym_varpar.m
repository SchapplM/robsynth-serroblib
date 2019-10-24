% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PPRPR3
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
%   Wie in S5PPRPR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:18
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PPRPR3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRPR3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:51
	% EndTime: 2019-10-24 10:18:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:51
	% EndTime: 2019-10-24 10:18:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:51
	% EndTime: 2019-10-24 10:18:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:51
	% EndTime: 2019-10-24 10:18:51
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (46->7), mult. (159->21), div. (18->4), fcn. (175->5), ass. (0->15)
	t43 = sin(pkin(7));
	t46 = sin(qJ(3));
	t47 = cos(qJ(3));
	t50 = cos(pkin(7)) * cos(pkin(8));
	t41 = t43 * t46 + t47 * t50;
	t38 = 0.1e1 / t41 ^ 2;
	t54 = qJD(3) * t38;
	t40 = -t43 * t47 + t46 * t50;
	t37 = t40 ^ 2;
	t34 = t37 * t38 + 0.1e1;
	t51 = t41 * t54;
	t52 = t40 / t41 * t54;
	t53 = (t37 * t52 + t40 * t51) / t34 ^ 2;
	t32 = 0.1e1 / t34;
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, -0.2e1 * t53 + 0.2e1 * (t32 * t51 + (t32 * t52 - t38 * t53) * t40) * t40, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:51
	% EndTime: 2019-10-24 10:18:51
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (116->8), mult. (159->21), div. (18->4), fcn. (175->5), ass. (0->16)
	t49 = qJ(3) + pkin(9);
	t47 = sin(t49);
	t48 = cos(t49);
	t50 = sin(pkin(7));
	t55 = cos(pkin(7)) * cos(pkin(8));
	t45 = t50 * t47 + t48 * t55;
	t42 = 0.1e1 / t45 ^ 2;
	t59 = qJD(3) * t42;
	t44 = t47 * t55 - t50 * t48;
	t41 = t44 ^ 2;
	t38 = t41 * t42 + 0.1e1;
	t56 = t45 * t59;
	t57 = t44 / t45 * t59;
	t58 = (t41 * t57 + t44 * t56) / t38 ^ 2;
	t36 = 0.1e1 / t38;
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, -0.2e1 * t58 + 0.2e1 * (t36 * t56 + (t36 * t57 - t42 * t58) * t44) * t44, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:51
	% EndTime: 2019-10-24 10:18:51
	% DurationCPUTime: 0.42s
	% Computational Cost: add. (1741->55), mult. (2271->130), div. (423->14), fcn. (2956->11), ass. (0->67)
	t145 = qJ(3) + pkin(9);
	t141 = sin(t145);
	t147 = sin(pkin(7));
	t148 = cos(pkin(8));
	t142 = cos(t145);
	t149 = cos(pkin(7));
	t169 = t142 * t149;
	t133 = t141 * t147 + t148 * t169;
	t150 = sin(qJ(5));
	t151 = cos(qJ(5));
	t146 = sin(pkin(8));
	t167 = t146 * t149;
	t158 = -t133 * t150 + t151 * t167;
	t178 = t158 * qJD(5);
	t165 = t149 * t141;
	t166 = t147 * t148;
	t131 = t142 * t166 - t165;
	t129 = t141 * t166 + t169;
	t168 = t146 * t141;
	t121 = atan2(-t129, t168);
	t117 = sin(t121);
	t118 = cos(t121);
	t104 = -t117 * t129 + t118 * t168;
	t101 = 0.1e1 / t104;
	t116 = t133 * t151 + t150 * t167;
	t112 = 0.1e1 / t116;
	t138 = 0.1e1 / t141;
	t102 = 0.1e1 / t104 ^ 2;
	t113 = 0.1e1 / t116 ^ 2;
	t139 = 0.1e1 / t141 ^ 2;
	t124 = t131 * qJD(3);
	t164 = qJD(3) * t142;
	t170 = t139 * t142;
	t162 = t129 * t170;
	t127 = t129 ^ 2;
	t144 = 0.1e1 / t146 ^ 2;
	t122 = t127 * t139 * t144 + 0.1e1;
	t119 = 0.1e1 / t122;
	t143 = 0.1e1 / t146;
	t171 = t119 * t143;
	t96 = (qJD(3) * t162 - t124 * t138) * t171;
	t93 = (-t129 * t96 + t146 * t164) * t118 + (-t96 * t168 - t124) * t117;
	t177 = t101 * t102 * t93;
	t111 = t158 ^ 2;
	t108 = t111 * t113 + 0.1e1;
	t132 = -t147 * t142 + t148 * t165;
	t125 = t132 * qJD(3);
	t109 = t116 * qJD(5) - t125 * t150;
	t172 = t113 * t158;
	t110 = -t125 * t151 + t178;
	t173 = t110 * t112 * t113;
	t176 = 0.1e1 / t108 ^ 2 * (-t109 * t172 - t111 * t173);
	t159 = -t131 * t138 + t162;
	t97 = t159 * t171;
	t175 = t129 * t97;
	t174 = t102 * t132;
	t163 = -0.2e1 * t176;
	t160 = -t112 * t150 - t151 * t172;
	t140 = t138 * t139;
	t128 = t132 ^ 2;
	t126 = t133 * qJD(3);
	t123 = t129 * qJD(3);
	t106 = 0.1e1 / t108;
	t100 = t102 * t128 + 0.1e1;
	t94 = (t142 * t146 - t175) * t118 + (-t97 * t168 - t131) * t117;
	t92 = (-0.2e1 * t159 / t122 ^ 2 * (t124 * t129 * t139 - t127 * t140 * t164) * t144 + (t124 * t170 + t123 * t138 + (t131 * t170 + (-0.2e1 * t140 * t142 ^ 2 - t138) * t129) * qJD(3)) * t119) * t143;
	t1 = [0, 0, t92, 0, 0; 0, 0, 0.2e1 * (-t101 * t133 + t94 * t174) / t100 ^ 2 * (t126 * t174 - t128 * t177) + (-t125 * t101 + (-t94 * t126 - t133 * t93) * t102 + (0.2e1 * t94 * t177 + (-(-t124 * t97 - t129 * t92 - t131 * t96 + (-t96 * t97 - qJD(3)) * t168) * t118 - (t96 * t175 + t123 + (-t141 * t92 + (-qJD(3) * t97 - t96) * t142) * t146) * t117) * t102) * t132) / t100, 0, 0; 0, 0, t160 * t132 * t163 + (t160 * t126 + ((-qJD(5) * t112 + 0.2e1 * t158 * t173) * t151 + (t109 * t151 + (t110 + t178) * t150) * t113) * t132) * t106, 0, t163 - 0.2e1 * (t106 * t109 * t113 - (-t106 * t173 - t113 * t176) * t158) * t158;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end