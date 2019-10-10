% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:12
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PPRRPR2_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRPR2_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_jacobigD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:34
	% EndTime: 2019-10-09 21:12:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:34
	% EndTime: 2019-10-09 21:12:34
	% DurationCPUTime: 0.05s
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
	% StartTime: 2019-10-09 21:12:34
	% EndTime: 2019-10-09 21:12:34
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (10->10), mult. (40->28), div. (0->0), fcn. (44->10), ass. (0->14)
	t156 = sin(pkin(11));
	t162 = cos(pkin(6));
	t167 = t156 * t162;
	t157 = sin(pkin(7));
	t158 = sin(pkin(6));
	t166 = t158 * t157;
	t160 = cos(pkin(11));
	t165 = t160 * t162;
	t164 = cos(qJ(3));
	t163 = sin(qJ(3));
	t161 = cos(pkin(7));
	t159 = cos(pkin(12));
	t155 = sin(pkin(12));
	t1 = [0, 0, 0, ((-t155 * t167 + t160 * t159) * t164 + ((-t160 * t155 - t159 * t167) * t161 + t156 * t166) * t163) * qJD(3), 0, 0; 0, 0, 0, ((t155 * t165 + t156 * t159) * t164 + ((-t156 * t155 + t159 * t165) * t161 - t160 * t166) * t163) * qJD(3), 0, 0; 0, 0, 0, (t157 * t162 * t163 + (t159 * t161 * t163 + t155 * t164) * t158) * qJD(3), 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:34
	% EndTime: 2019-10-09 21:12:34
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (41->23), mult. (141->56), div. (0->0), fcn. (164->12), ass. (0->29)
	t215 = sin(pkin(11));
	t221 = cos(pkin(6));
	t237 = t215 * t221;
	t216 = sin(pkin(7));
	t217 = sin(pkin(6));
	t236 = t216 * t217;
	t235 = t216 * t221;
	t220 = cos(pkin(7));
	t234 = t217 * t220;
	t218 = cos(pkin(12));
	t233 = t218 * t220;
	t219 = cos(pkin(11));
	t232 = t219 * t221;
	t224 = cos(qJ(4));
	t231 = qJD(3) * t224;
	t214 = sin(pkin(12));
	t210 = -t215 * t214 + t218 * t232;
	t230 = t210 * t220 - t219 * t236;
	t212 = -t219 * t214 - t218 * t237;
	t229 = t212 * t220 + t215 * t236;
	t211 = t214 * t232 + t215 * t218;
	t223 = sin(qJ(3));
	t225 = cos(qJ(3));
	t228 = t211 * t225 + t230 * t223;
	t213 = -t214 * t237 + t219 * t218;
	t227 = t213 * t225 + t229 * t223;
	t226 = t223 * t235 + (t214 * t225 + t223 * t233) * t217;
	t222 = sin(qJ(4));
	t1 = [0, 0, 0, t227 * qJD(3), 0, (-t227 * t222 + (-t212 * t216 + t215 * t234) * t224) * qJD(4) + (-t213 * t223 + t229 * t225) * t231; 0, 0, 0, t228 * qJD(3), 0, (-t228 * t222 + (-t210 * t216 - t219 * t234) * t224) * qJD(4) + (-t211 * t223 + t230 * t225) * t231; 0, 0, 0, t226 * qJD(3), 0, (-t226 * t222 + (-t218 * t236 + t221 * t220) * t224) * qJD(4) + (t225 * t235 + (-t214 * t223 + t225 * t233) * t217) * t231;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end