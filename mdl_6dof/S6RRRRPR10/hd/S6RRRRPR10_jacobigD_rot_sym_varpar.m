% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR10
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:46
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRPR10_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR10_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
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
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (12->6), mult. (48->16), div. (0->0), fcn. (48->6), ass. (0->15)
	t128 = sin(qJ(2));
	t129 = sin(qJ(1));
	t136 = t128 * t129;
	t131 = cos(qJ(1));
	t135 = t128 * t131;
	t130 = cos(qJ(2));
	t134 = t129 * t130;
	t133 = t130 * t131;
	t126 = sin(pkin(6));
	t132 = qJD(1) * t126;
	t127 = cos(pkin(6));
	t125 = t126 * qJD(2) * t128;
	t124 = (t127 * t135 + t134) * qJD(2) + (t127 * t134 + t135) * qJD(1);
	t123 = (-t127 * t136 + t133) * qJD(2) + (t127 * t133 - t136) * qJD(1);
	t1 = [0, t131 * t132, t123, t123, 0, 0; 0, t129 * t132, t124, t124, 0, 0; 0, 0, t125, t125, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:06
	% EndTime: 2019-10-10 12:46:06
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (12->6), mult. (48->16), div. (0->0), fcn. (48->6), ass. (0->15)
	t149 = sin(qJ(2));
	t150 = sin(qJ(1));
	t157 = t149 * t150;
	t152 = cos(qJ(1));
	t156 = t149 * t152;
	t151 = cos(qJ(2));
	t155 = t150 * t151;
	t154 = t151 * t152;
	t147 = sin(pkin(6));
	t153 = qJD(1) * t147;
	t148 = cos(pkin(6));
	t146 = t147 * qJD(2) * t149;
	t145 = (t148 * t156 + t155) * qJD(2) + (t148 * t155 + t156) * qJD(1);
	t144 = (-t148 * t157 + t154) * qJD(2) + (t148 * t154 - t157) * qJD(1);
	t1 = [0, t152 * t153, t144, t144, 0, 0; 0, t150 * t153, t145, t145, 0, 0; 0, 0, t146, t146, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:06
	% EndTime: 2019-10-10 12:46:06
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (45->18), mult. (100->40), div. (0->0), fcn. (102->8), ass. (0->27)
	t196 = sin(pkin(6));
	t199 = sin(qJ(1));
	t214 = t196 * t199;
	t201 = cos(qJ(1));
	t213 = t196 * t201;
	t198 = sin(qJ(2));
	t212 = t198 * t199;
	t211 = t198 * t201;
	t200 = cos(qJ(2));
	t210 = t199 * t200;
	t209 = t201 * t200;
	t208 = qJD(1) * t196;
	t195 = qJ(3) + qJ(4);
	t193 = cos(t195);
	t207 = qJD(2) * t193;
	t206 = qJD(2) * t196;
	t197 = cos(pkin(6));
	t205 = t197 * t209 - t212;
	t204 = t197 * t210 + t211;
	t203 = t197 * t211 + t210;
	t202 = -t197 * t212 + t209;
	t194 = qJD(3) + qJD(4);
	t192 = sin(t195);
	t191 = t198 * t206;
	t190 = qJD(1) * t204 + qJD(2) * t203;
	t189 = qJD(1) * t205 + qJD(2) * t202;
	t1 = [0, t201 * t208, t189, t189, 0, (-t192 * t202 + t193 * t214) * t194 - t204 * t207 + (t192 * t213 - t193 * t203) * qJD(1); 0, t199 * t208, t190, t190, 0, (-t192 * t203 - t193 * t213) * t194 + t205 * t207 + (t192 * t214 + t193 * t202) * qJD(1); 0, 0, t191, t191, 0, -t196 * t198 * t194 * t192 + (t194 * t197 + t200 * t206) * t193;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end