% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:08
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PPRPRR1_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRPRR1_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobigD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:44
	% EndTime: 2019-10-09 21:08:44
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (19->12), mult. (73->29), div. (0->0), fcn. (77->12), ass. (0->18)
	t168 = sin(pkin(11));
	t175 = cos(pkin(6));
	t181 = t168 * t175;
	t166 = sin(pkin(13));
	t171 = cos(pkin(13));
	t176 = sin(qJ(3));
	t177 = cos(qJ(3));
	t178 = qJD(3) * (-t166 * t177 - t171 * t176);
	t163 = sin(pkin(7)) * t178;
	t170 = sin(pkin(6));
	t180 = t170 * t163;
	t173 = cos(pkin(11));
	t179 = t173 * t175;
	t172 = cos(pkin(12));
	t167 = sin(pkin(12));
	t165 = (t166 * t176 - t171 * t177) * qJD(3);
	t164 = cos(pkin(7)) * t178;
	t1 = [0, 0, 0, 0, -(-t167 * t181 + t173 * t172) * t165 - (-t173 * t167 - t172 * t181) * t164 - t168 * t180, 0; 0, 0, 0, 0, -(t167 * t179 + t168 * t172) * t165 - (-t168 * t167 + t172 * t179) * t164 + t173 * t180, 0; 0, 0, 0, 0, -t175 * t163 + (-t164 * t172 - t165 * t167) * t170, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:44
	% EndTime: 2019-10-09 21:08:44
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (68->33), mult. (231->74), div. (0->0), fcn. (263->14), ass. (0->34)
	t241 = sin(pkin(13));
	t246 = cos(pkin(13));
	t252 = sin(qJ(3));
	t254 = cos(qJ(3));
	t256 = t252 * t241 - t254 * t246;
	t263 = t256 * qJD(3);
	t243 = sin(pkin(11));
	t245 = sin(pkin(6));
	t262 = t243 * t245;
	t250 = cos(pkin(6));
	t261 = t243 * t250;
	t248 = cos(pkin(11));
	t260 = t245 * t248;
	t249 = cos(pkin(7));
	t259 = t245 * t249;
	t258 = t248 * t250;
	t257 = t241 * t254 + t246 * t252;
	t239 = t257 * qJD(3);
	t253 = cos(qJ(5));
	t251 = sin(qJ(5));
	t247 = cos(pkin(12));
	t244 = sin(pkin(7));
	t242 = sin(pkin(12));
	t237 = -t242 * t261 + t248 * t247;
	t236 = -t248 * t242 - t247 * t261;
	t235 = t242 * t258 + t243 * t247;
	t234 = -t243 * t242 + t247 * t258;
	t233 = t257 * t249;
	t232 = t257 * t244;
	t231 = t249 * t263;
	t230 = t249 * t239;
	t229 = t244 * t263;
	t228 = t244 * t239;
	t1 = [0, 0, 0, 0, t228 * t262 + t236 * t230 - t237 * t263, (-t229 * t262 - t236 * t231 - t237 * t239) * t251 + ((t232 * t262 + t236 * t233 - t237 * t256) * t253 + (-t236 * t244 + t243 * t259) * t251) * qJD(5); 0, 0, 0, 0, -t228 * t260 + t234 * t230 - t235 * t263, (t229 * t260 - t234 * t231 - t235 * t239) * t251 + ((-t232 * t260 + t234 * t233 - t235 * t256) * t253 + (-t234 * t244 - t248 * t259) * t251) * qJD(5); 0, 0, 0, 0, t250 * t228 + (t230 * t247 - t242 * t263) * t245, (-t250 * t229 + (-t231 * t247 - t239 * t242) * t245) * t251 + ((t250 * t232 + (t233 * t247 - t242 * t256) * t245) * t253 + (-t245 * t247 * t244 + t250 * t249) * t251) * qJD(5);];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end