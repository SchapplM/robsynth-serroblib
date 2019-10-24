% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PPPRR2
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
%   Wie in S5PPPRR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:17
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PPPRR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPPRR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:41
	% EndTime: 2019-10-24 10:17:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:41
	% EndTime: 2019-10-24 10:17:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:41
	% EndTime: 2019-10-24 10:17:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:41
	% EndTime: 2019-10-24 10:17:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:41
	% EndTime: 2019-10-24 10:17:41
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (81->8), mult. (264->24), div. (18->4), fcn. (315->8), ass. (0->16)
	t57 = cos(pkin(7));
	t54 = t57 * cos(pkin(8)) * cos(pkin(9)) + sin(pkin(7)) * sin(pkin(9));
	t58 = sin(qJ(4));
	t59 = cos(qJ(4));
	t63 = sin(pkin(8)) * t57;
	t53 = t54 * t59 + t58 * t63;
	t50 = 0.1e1 / t53 ^ 2;
	t67 = qJD(4) * t50;
	t62 = -t54 * t58 + t59 * t63;
	t49 = t62 ^ 2;
	t46 = t49 * t50 + 0.1e1;
	t64 = t53 * t67;
	t65 = t62 / t53 * t67;
	t66 = (-t49 * t65 - t62 * t64) / t46 ^ 2;
	t44 = 0.1e1 / t46;
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, -0.2e1 * t66 - 0.2e1 * (t44 * t64 - (-t44 * t65 - t50 * t66) * t62) * t62, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:41
	% EndTime: 2019-10-24 10:17:41
	% DurationCPUTime: 0.43s
	% Computational Cost: add. (1445->58), mult. (4330->131), div. (281->12), fcn. (5680->13), ass. (0->71)
	t175 = cos(pkin(8));
	t174 = cos(pkin(9));
	t205 = sin(pkin(7));
	t191 = t205 * t174;
	t173 = sin(pkin(9));
	t206 = cos(pkin(7));
	t192 = t206 * t173;
	t163 = t175 * t191 - t192;
	t177 = sin(qJ(4));
	t204 = sin(pkin(8));
	t189 = t205 * t204;
	t207 = cos(qJ(4));
	t155 = t163 * t177 - t207 * t189;
	t193 = t174 * t204;
	t169 = t175 * t207 + t177 * t193;
	t148 = atan2(-t155, t169);
	t143 = sin(t148);
	t144 = cos(t148);
	t130 = -t143 * t155 + t144 * t169;
	t127 = 0.1e1 / t130;
	t165 = t206 * t175 * t174 + t205 * t173;
	t190 = t206 * t204;
	t159 = t165 * t207 + t177 * t190;
	t164 = t175 * t192 - t191;
	t176 = sin(qJ(5));
	t178 = cos(qJ(5));
	t142 = t159 * t178 + t164 * t176;
	t138 = 0.1e1 / t142;
	t166 = 0.1e1 / t169;
	t128 = 0.1e1 / t130 ^ 2;
	t139 = 0.1e1 / t142 ^ 2;
	t167 = 0.1e1 / t169 ^ 2;
	t153 = t155 ^ 2;
	t147 = t153 * t167 + 0.1e1;
	t145 = 0.1e1 / t147;
	t157 = t163 * t207 + t177 * t189;
	t150 = t157 * qJD(4);
	t170 = -t175 * t177 + t207 * t193;
	t162 = t170 * qJD(4);
	t198 = t155 * t167;
	t122 = (-t150 * t166 + t162 * t198) * t145;
	t188 = -t143 * t169 - t144 * t155;
	t119 = t188 * t122 - t143 * t150 + t144 * t162;
	t203 = t119 * t127 * t128;
	t141 = t159 * t176 - t164 * t178;
	t137 = t141 ^ 2;
	t134 = t137 * t139 + 0.1e1;
	t185 = -t165 * t177 + t207 * t190;
	t151 = t185 * qJD(4);
	t135 = t142 * qJD(5) + t151 * t176;
	t199 = t139 * t141;
	t195 = qJD(5) * t141;
	t136 = t151 * t178 - t195;
	t200 = t136 * t138 * t139;
	t202 = (t135 * t199 - t137 * t200) / t134 ^ 2;
	t201 = t128 * t185;
	t197 = t155 * t170;
	t196 = t162 * t166 * t167;
	t194 = -0.2e1 * t202;
	t187 = -t138 * t176 + t178 * t199;
	t186 = -t157 * t166 + t167 * t197;
	t161 = t169 * qJD(4);
	t154 = t185 ^ 2;
	t152 = t159 * qJD(4);
	t149 = t155 * qJD(4);
	t132 = 0.1e1 / t134;
	t126 = t154 * t128 + 0.1e1;
	t123 = t186 * t145;
	t120 = t188 * t123 - t143 * t157 + t144 * t170;
	t118 = -0.2e1 * t186 / t147 ^ 2 * (t150 * t198 - t153 * t196) + (-0.2e1 * t196 * t197 + t149 * t166 + (t150 * t170 - t155 * t161 + t157 * t162) * t167) * t145;
	t1 = [0, 0, 0, t118, 0; 0, 0, 0, 0.2e1 * (-t120 * t201 - t127 * t159) / t126 ^ 2 * (-t152 * t201 - t154 * t203) + (t151 * t127 + (-t159 * t119 - t120 * t152) * t128 - (0.2e1 * t120 * t203 + (-(-t118 * t155 - t123 * t150 - t161 + (-t123 * t169 - t157) * t122) * t144 - (-t118 * t169 - t123 * t162 + t149 + (t123 * t155 - t170) * t122) * t143) * t128) * t185) / t126, 0; 0, 0, 0, -t187 * t185 * t194 + (t187 * t152 - ((-qJD(5) * t138 - 0.2e1 * t141 * t200) * t178 + (t135 * t178 + (t136 - t195) * t176) * t139) * t185) * t132, t194 + 0.2e1 * (t132 * t135 * t139 + (-t132 * t200 - t139 * t202) * t141) * t141;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end