% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR12
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:12
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRR12_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR12_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR12_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_jacobigD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:12:18
	% EndTime: 2019-10-10 12:12:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:12:18
	% EndTime: 2019-10-10 12:12:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:12:18
	% EndTime: 2019-10-10 12:12:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:12:19
	% EndTime: 2019-10-10 12:12:19
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
	t103 = sin(qJ(2));
	t104 = sin(qJ(1));
	t111 = t103 * t104;
	t106 = cos(qJ(1));
	t110 = t103 * t106;
	t105 = cos(qJ(2));
	t109 = t104 * t105;
	t108 = t105 * t106;
	t101 = sin(pkin(6));
	t107 = qJD(1) * t101;
	t102 = cos(pkin(6));
	t1 = [0, t106 * t107, (-t102 * t111 + t108) * qJD(2) + (t102 * t108 - t111) * qJD(1), 0, 0, 0; 0, t104 * t107, (t102 * t110 + t109) * qJD(2) + (t102 * t109 + t110) * qJD(1), 0, 0, 0; 0, 0, t101 * qJD(2) * t103, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:12:19
	% EndTime: 2019-10-10 12:12:19
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
	t133 = sin(qJ(2));
	t134 = sin(qJ(1));
	t141 = t133 * t134;
	t136 = cos(qJ(1));
	t140 = t133 * t136;
	t135 = cos(qJ(2));
	t139 = t134 * t135;
	t138 = t135 * t136;
	t131 = sin(pkin(6));
	t137 = qJD(1) * t131;
	t132 = cos(pkin(6));
	t1 = [0, t136 * t137, (-t132 * t141 + t138) * qJD(2) + (t132 * t138 - t141) * qJD(1), 0, 0, 0; 0, t134 * t137, (t132 * t140 + t139) * qJD(2) + (t132 * t139 + t140) * qJD(1), 0, 0, 0; 0, 0, t131 * qJD(2) * t133, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:12:19
	% EndTime: 2019-10-10 12:12:19
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (22->16), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->21)
	t158 = sin(pkin(6));
	t163 = cos(qJ(3));
	t177 = t158 * t163;
	t165 = cos(qJ(1));
	t176 = t158 * t165;
	t161 = sin(qJ(2));
	t162 = sin(qJ(1));
	t175 = t161 * t162;
	t174 = t161 * t165;
	t164 = cos(qJ(2));
	t173 = t162 * t164;
	t172 = t165 * t164;
	t171 = qJD(1) * t158;
	t160 = sin(qJ(3));
	t170 = qJD(2) * t160;
	t159 = cos(pkin(6));
	t169 = t159 * t172 - t175;
	t168 = t159 * t173 + t174;
	t167 = t159 * t174 + t173;
	t166 = -t159 * t175 + t172;
	t1 = [0, t165 * t171, t169 * qJD(1) + t166 * qJD(2), 0, (t162 * t158 * t160 + t166 * t163) * qJD(3) - t168 * t170 + (-t167 * t160 - t163 * t176) * qJD(1), 0; 0, t162 * t171, t168 * qJD(1) + t167 * qJD(2), 0, (-t160 * t176 + t167 * t163) * qJD(3) + t169 * t170 + (t166 * t160 - t162 * t177) * qJD(1), 0; 0, 0, t158 * qJD(2) * t161, 0, t158 * t164 * t170 + (t159 * t160 + t161 * t177) * qJD(3), 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:12:19
	% EndTime: 2019-10-10 12:12:19
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (38->16), mult. (130->40), div. (0->0), fcn. (134->8), ass. (0->24)
	t179 = sin(pkin(6));
	t184 = cos(qJ(3));
	t198 = t179 * t184;
	t186 = cos(qJ(1));
	t197 = t179 * t186;
	t182 = sin(qJ(2));
	t183 = sin(qJ(1));
	t196 = t182 * t183;
	t195 = t182 * t186;
	t185 = cos(qJ(2));
	t194 = t183 * t185;
	t193 = t186 * t185;
	t192 = qJD(1) * t179;
	t181 = sin(qJ(3));
	t191 = qJD(2) * t181;
	t180 = cos(pkin(6));
	t190 = t180 * t193 - t196;
	t189 = t180 * t194 + t195;
	t188 = t180 * t195 + t194;
	t187 = -t180 * t196 + t193;
	t178 = t179 * t185 * t191 + (t180 * t181 + t182 * t198) * qJD(3);
	t177 = (-t181 * t197 + t188 * t184) * qJD(3) + t190 * t191 + (t187 * t181 - t183 * t198) * qJD(1);
	t176 = (t183 * t179 * t181 + t187 * t184) * qJD(3) - t189 * t191 + (-t188 * t181 - t184 * t197) * qJD(1);
	t1 = [0, t186 * t192, t190 * qJD(1) + t187 * qJD(2), 0, t176, t176; 0, t183 * t192, t189 * qJD(1) + t188 * qJD(2), 0, t177, t177; 0, 0, t179 * qJD(2) * t182, 0, t178, t178;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end