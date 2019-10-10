% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PPRRPR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:10
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PPRRPR1_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRPR1_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_jacobigD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (10->10), mult. (40->28), div. (0->0), fcn. (44->10), ass. (0->14)
	t136 = sin(pkin(11));
	t142 = cos(pkin(6));
	t147 = t136 * t142;
	t137 = sin(pkin(7));
	t138 = sin(pkin(6));
	t146 = t138 * t137;
	t140 = cos(pkin(11));
	t145 = t140 * t142;
	t144 = cos(qJ(3));
	t143 = sin(qJ(3));
	t141 = cos(pkin(7));
	t139 = cos(pkin(12));
	t135 = sin(pkin(12));
	t1 = [0, 0, 0, ((-t135 * t147 + t140 * t139) * t144 + ((-t140 * t135 - t139 * t147) * t141 + t136 * t146) * t143) * qJD(3), 0, 0; 0, 0, 0, ((t135 * t145 + t136 * t139) * t144 + ((-t136 * t135 + t139 * t145) * t141 - t140 * t146) * t143) * qJD(3), 0, 0; 0, 0, 0, (t137 * t142 * t143 + (t139 * t141 * t143 + t135 * t144) * t138) * qJD(3), 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (10->10), mult. (40->28), div. (0->0), fcn. (44->10), ass. (0->14)
	t177 = sin(pkin(11));
	t183 = cos(pkin(6));
	t188 = t177 * t183;
	t178 = sin(pkin(7));
	t179 = sin(pkin(6));
	t187 = t179 * t178;
	t181 = cos(pkin(11));
	t186 = t181 * t183;
	t185 = cos(qJ(3));
	t184 = sin(qJ(3));
	t182 = cos(pkin(7));
	t180 = cos(pkin(12));
	t176 = sin(pkin(12));
	t1 = [0, 0, 0, ((-t176 * t188 + t181 * t180) * t185 + ((-t181 * t176 - t180 * t188) * t182 + t177 * t187) * t184) * qJD(3), 0, 0; 0, 0, 0, ((t176 * t186 + t177 * t180) * t185 + ((-t177 * t176 + t180 * t186) * t182 - t181 * t187) * t184) * qJD(3), 0, 0; 0, 0, 0, (t178 * t183 * t184 + (t180 * t182 * t184 + t176 * t185) * t179) * qJD(3), 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:41
	% EndTime: 2019-10-09 21:10:41
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (41->23), mult. (141->56), div. (0->0), fcn. (164->12), ass. (0->29)
	t223 = sin(pkin(11));
	t229 = cos(pkin(6));
	t245 = t223 * t229;
	t224 = sin(pkin(7));
	t225 = sin(pkin(6));
	t244 = t224 * t225;
	t243 = t224 * t229;
	t228 = cos(pkin(7));
	t242 = t225 * t228;
	t226 = cos(pkin(12));
	t241 = t226 * t228;
	t227 = cos(pkin(11));
	t240 = t227 * t229;
	t230 = sin(qJ(4));
	t239 = qJD(3) * t230;
	t222 = sin(pkin(12));
	t218 = -t223 * t222 + t226 * t240;
	t238 = t218 * t228 - t227 * t244;
	t220 = -t227 * t222 - t226 * t245;
	t237 = t220 * t228 + t223 * t244;
	t219 = t222 * t240 + t223 * t226;
	t231 = sin(qJ(3));
	t233 = cos(qJ(3));
	t236 = t219 * t233 + t238 * t231;
	t221 = -t222 * t245 + t227 * t226;
	t235 = t221 * t233 + t237 * t231;
	t234 = t231 * t243 + (t222 * t233 + t231 * t241) * t225;
	t232 = cos(qJ(4));
	t1 = [0, 0, 0, t235 * qJD(3), 0, (t235 * t232 + (-t220 * t224 + t223 * t242) * t230) * qJD(4) + (-t221 * t231 + t237 * t233) * t239; 0, 0, 0, t236 * qJD(3), 0, (t236 * t232 + (-t218 * t224 - t227 * t242) * t230) * qJD(4) + (-t219 * t231 + t238 * t233) * t239; 0, 0, 0, t234 * qJD(3), 0, (t234 * t232 + (-t226 * t244 + t229 * t228) * t230) * qJD(4) + (t233 * t243 + (-t222 * t231 + t233 * t241) * t225) * t239;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end