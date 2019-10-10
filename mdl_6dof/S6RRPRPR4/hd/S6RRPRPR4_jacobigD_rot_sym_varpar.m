% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:09
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRPR4_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR4_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_jacobigD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:42
	% EndTime: 2019-10-10 10:09:42
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t79 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t79, 0, 0, 0, 0; 0, sin(qJ(1)) * t79, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:42
	% EndTime: 2019-10-10 10:09:42
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (15->8), mult. (54->22), div. (0->0), fcn. (58->8), ass. (0->16)
	t133 = sin(pkin(11));
	t135 = cos(pkin(11));
	t137 = sin(qJ(2));
	t139 = cos(qJ(2));
	t142 = t139 * t133 + t137 * t135;
	t144 = qJD(2) * t142;
	t134 = sin(pkin(6));
	t143 = qJD(1) * t134;
	t141 = t133 * t137 - t135 * t139;
	t140 = cos(qJ(1));
	t138 = sin(qJ(1));
	t136 = cos(pkin(6));
	t131 = t141 * qJD(2);
	t130 = t141 * t136;
	t129 = t136 * t144;
	t1 = [0, t140 * t143, 0, -t138 * t129 - t140 * t131 + (-t130 * t140 - t138 * t142) * qJD(1), 0, 0; 0, t138 * t143, 0, t140 * t129 - t138 * t131 + (-t130 * t138 + t140 * t142) * qJD(1), 0, 0; 0, 0, 0, t134 * t144, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:42
	% EndTime: 2019-10-10 10:09:42
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (15->8), mult. (54->22), div. (0->0), fcn. (58->8), ass. (0->16)
	t143 = sin(pkin(11));
	t145 = cos(pkin(11));
	t147 = sin(qJ(2));
	t149 = cos(qJ(2));
	t152 = t149 * t143 + t147 * t145;
	t154 = qJD(2) * t152;
	t144 = sin(pkin(6));
	t153 = qJD(1) * t144;
	t151 = t143 * t147 - t145 * t149;
	t150 = cos(qJ(1));
	t148 = sin(qJ(1));
	t146 = cos(pkin(6));
	t141 = t151 * qJD(2);
	t140 = t151 * t146;
	t139 = t146 * t154;
	t1 = [0, t150 * t153, 0, -t148 * t139 - t150 * t141 + (-t140 * t150 - t148 * t152) * qJD(1), 0, 0; 0, t148 * t153, 0, t150 * t139 - t148 * t141 + (-t140 * t148 + t150 * t152) * qJD(1), 0, 0; 0, 0, 0, t144 * t154, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:42
	% EndTime: 2019-10-10 10:09:42
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (56->23), mult. (147->54), div. (0->0), fcn. (162->10), ass. (0->25)
	t213 = sin(pkin(11));
	t215 = cos(pkin(11));
	t217 = sin(qJ(2));
	t219 = cos(qJ(2));
	t222 = t217 * t213 - t219 * t215;
	t229 = t222 * qJD(2);
	t223 = t219 * t213 + t217 * t215;
	t207 = t223 * qJD(2);
	t214 = sin(pkin(6));
	t218 = sin(qJ(1));
	t228 = t214 * t218;
	t220 = cos(qJ(1));
	t227 = t214 * t220;
	t226 = qJD(1) * t214;
	t216 = cos(pkin(6));
	t205 = t223 * t216;
	t225 = t220 * t205 - t218 * t222;
	t224 = -t218 * t205 - t220 * t222;
	t212 = qJ(4) + pkin(12);
	t211 = cos(t212);
	t210 = sin(t212);
	t204 = t222 * t216;
	t203 = t216 * t229;
	t202 = t216 * t207;
	t1 = [0, t220 * t226, 0, -t218 * t202 - t220 * t229 + (-t204 * t220 - t218 * t223) * qJD(1), 0, (t218 * t203 - t220 * t207) * t210 + (t210 * t228 + t224 * t211) * qJD(4) + (-t225 * t210 - t211 * t227) * qJD(1); 0, t218 * t226, 0, t220 * t202 - t218 * t229 + (-t204 * t218 + t220 * t223) * qJD(1), 0, (-t220 * t203 - t218 * t207) * t210 + (-t210 * t227 + t225 * t211) * qJD(4) + (t224 * t210 - t211 * t228) * qJD(1); 0, 0, 0, t214 * t207, 0, t216 * qJD(4) * t210 + (t223 * qJD(4) * t211 - t210 * t229) * t214;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end