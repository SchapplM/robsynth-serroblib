% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:46
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPPRPR8_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR8_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:46:04
	% EndTime: 2019-10-09 23:46:04
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:46:04
	% EndTime: 2019-10-09 23:46:05
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:46:05
	% EndTime: 2019-10-09 23:46:05
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t13 = qJD(1) * sin(qJ(1));
	t11 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; t13, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t13, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:46:04
	% EndTime: 2019-10-09 23:46:04
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (4->4), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t24 = qJD(1) * sin(qJ(1));
	t23 = qJD(1) * cos(qJ(1));
	t20 = cos(pkin(9));
	t19 = sin(pkin(9));
	t1 = [-t19 * t24, 0, 0, 0, 0, 0; t19 * t23, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t20 * t24, 0, 0, 0, 0, 0; t20 * t23, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t23, 0, 0, 0, 0, 0; -t24, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:46:05
	% EndTime: 2019-10-09 23:46:05
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t44 = sin(qJ(1));
	t49 = qJD(1) * t44;
	t45 = cos(qJ(1));
	t48 = qJD(1) * t45;
	t47 = qJD(4) * t44;
	t46 = qJD(4) * t45;
	t43 = pkin(9) + qJ(4);
	t42 = cos(t43);
	t41 = sin(t43);
	t40 = -t41 * t47 + t42 * t48;
	t39 = t41 * t48 + t42 * t47;
	t38 = t41 * t46 + t42 * t49;
	t37 = -t41 * t49 + t42 * t46;
	t1 = [t37, 0, 0, t40, 0, 0; t39, 0, 0, t38, 0, 0; 0, 0, 0, -qJD(4) * t42, 0, 0; -t38, 0, 0, -t39, 0, 0; t40, 0, 0, t37, 0, 0; 0, 0, 0, qJD(4) * t41, 0, 0; -t48, 0, 0, 0, 0, 0; -t49, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:46:05
	% EndTime: 2019-10-09 23:46:05
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t172 = sin(qJ(1));
	t177 = qJD(1) * t172;
	t173 = cos(qJ(1));
	t176 = qJD(1) * t173;
	t175 = qJD(4) * t172;
	t174 = qJD(4) * t173;
	t171 = pkin(9) + qJ(4);
	t170 = cos(t171);
	t169 = sin(t171);
	t168 = t169 * t175 - t170 * t176;
	t167 = t169 * t176 + t170 * t175;
	t166 = t169 * t174 + t170 * t177;
	t165 = t169 * t177 - t170 * t174;
	t1 = [-t176, 0, 0, 0, 0, 0; -t177, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t165, 0, 0, t168, 0, 0; -t167, 0, 0, -t166, 0, 0; 0, 0, 0, qJD(4) * t170, 0, 0; t166, 0, 0, t167, 0, 0; t168, 0, 0, t165, 0, 0; 0, 0, 0, -qJD(4) * t169, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:46:06
	% EndTime: 2019-10-09 23:46:06
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (102->29), mult. (173->57), div. (0->0), fcn. (173->6), ass. (0->31)
	t249 = cos(qJ(6));
	t250 = cos(qJ(1));
	t269 = t249 * t250;
	t248 = sin(qJ(1));
	t268 = qJD(1) * t248;
	t267 = qJD(1) * t250;
	t247 = sin(qJ(6));
	t266 = qJD(4) * t247;
	t265 = qJD(4) * t248;
	t264 = qJD(4) * t249;
	t263 = qJD(4) * t250;
	t262 = qJD(6) * t247;
	t261 = qJD(6) * t249;
	t260 = qJD(6) * t250;
	t246 = pkin(9) + qJ(4);
	t244 = sin(t246);
	t259 = t244 * t264;
	t258 = t244 * t265;
	t245 = cos(t246);
	t257 = t245 * t265;
	t256 = t244 * t263;
	t255 = t245 * t263;
	t254 = qJD(6) * t245 + qJD(1);
	t253 = qJD(1) * t245 + qJD(6);
	t252 = t254 * t247;
	t251 = t253 * t248 + t256;
	t243 = t251 * t247 - t254 * t269;
	t242 = t251 * t249 + t250 * t252;
	t241 = t254 * t249 * t248 + (t253 * t250 - t258) * t247;
	t240 = -t253 * t269 + (t252 + t259) * t248;
	t1 = [t243, 0, 0, t247 * t257 + (t247 * t267 + t248 * t261) * t244, 0, t240; -t241, 0, 0, -t247 * t255 + (t247 * t268 - t249 * t260) * t244, 0, -t242; 0, 0, 0, -t244 * t266 + t245 * t261, 0, -t244 * t262 + t245 * t264; t242, 0, 0, t249 * t257 + (-t248 * t262 + t249 * t267) * t244, 0, t241; t240, 0, 0, -t249 * t255 + (t247 * t260 + t249 * t268) * t244, 0, t243; 0, 0, 0, -t245 * t262 - t259, 0, -t244 * t261 - t245 * t266; -t244 * t268 + t255, 0, 0, t245 * t267 - t258, 0, 0; t244 * t267 + t257, 0, 0, t245 * t268 + t256, 0, 0; 0, 0, 0, -qJD(4) * t245, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end