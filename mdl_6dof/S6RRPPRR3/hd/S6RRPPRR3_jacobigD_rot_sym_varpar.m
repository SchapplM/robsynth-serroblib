% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:39
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPPRR3_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR3_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobigD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:19
	% EndTime: 2019-10-10 09:39:19
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t79 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t79, 0, 0, 0, 0; 0, sin(qJ(1)) * t79, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:19
	% EndTime: 2019-10-10 09:39:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t113 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t113, 0, 0, 0, 0; 0, sin(qJ(1)) * t113, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:19
	% EndTime: 2019-10-10 09:39:19
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (15->8), mult. (54->22), div. (0->0), fcn. (58->8), ass. (0->16)
	t142 = sin(pkin(11));
	t144 = cos(pkin(11));
	t146 = sin(qJ(2));
	t148 = cos(qJ(2));
	t151 = t148 * t142 + t146 * t144;
	t153 = qJD(2) * t151;
	t143 = sin(pkin(6));
	t152 = qJD(1) * t143;
	t150 = t142 * t146 - t144 * t148;
	t149 = cos(qJ(1));
	t147 = sin(qJ(1));
	t145 = cos(pkin(6));
	t140 = t150 * qJD(2);
	t139 = t150 * t145;
	t138 = t145 * t153;
	t1 = [0, t149 * t152, 0, 0, -t147 * t138 - t149 * t140 + (-t139 * t149 - t147 * t151) * qJD(1), 0; 0, t147 * t152, 0, 0, t149 * t138 - t147 * t140 + (-t139 * t147 + t149 * t151) * qJD(1), 0; 0, 0, 0, 0, t143 * t153, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:19
	% EndTime: 2019-10-10 09:39:20
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (56->23), mult. (147->54), div. (0->0), fcn. (162->10), ass. (0->25)
	t212 = sin(pkin(11));
	t214 = cos(pkin(11));
	t216 = sin(qJ(2));
	t218 = cos(qJ(2));
	t221 = t216 * t212 - t218 * t214;
	t228 = t221 * qJD(2);
	t222 = t218 * t212 + t216 * t214;
	t206 = t222 * qJD(2);
	t213 = sin(pkin(6));
	t217 = sin(qJ(1));
	t227 = t213 * t217;
	t219 = cos(qJ(1));
	t226 = t213 * t219;
	t225 = qJD(1) * t213;
	t215 = cos(pkin(6));
	t204 = t222 * t215;
	t224 = t219 * t204 - t217 * t221;
	t223 = -t217 * t204 - t219 * t221;
	t211 = pkin(12) + qJ(5);
	t210 = cos(t211);
	t209 = sin(t211);
	t203 = t221 * t215;
	t202 = t215 * t228;
	t201 = t215 * t206;
	t1 = [0, t219 * t225, 0, 0, -t217 * t201 - t219 * t228 + (-t203 * t219 - t217 * t222) * qJD(1), (t217 * t202 - t219 * t206) * t209 + (t209 * t227 + t223 * t210) * qJD(5) + (-t224 * t209 - t210 * t226) * qJD(1); 0, t217 * t225, 0, 0, t219 * t201 - t217 * t228 + (-t203 * t217 + t219 * t222) * qJD(1), (-t219 * t202 - t217 * t206) * t209 + (-t209 * t226 + t224 * t210) * qJD(5) + (t223 * t209 - t210 * t227) * qJD(1); 0, 0, 0, 0, t213 * t206, t215 * qJD(5) * t209 + (t222 * qJD(5) * t210 - t209 * t228) * t213;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end