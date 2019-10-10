% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:25
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRPPR7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:25:42
	% EndTime: 2019-10-10 00:25:42
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:25:42
	% EndTime: 2019-10-10 00:25:42
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
	% StartTime: 2019-10-10 00:25:42
	% EndTime: 2019-10-10 00:25:42
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
	% StartTime: 2019-10-10 00:25:42
	% EndTime: 2019-10-10 00:25:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t35 = sin(qJ(1));
	t42 = qJD(1) * t35;
	t37 = cos(qJ(1));
	t41 = qJD(1) * t37;
	t34 = sin(qJ(3));
	t40 = qJD(3) * t34;
	t36 = cos(qJ(3));
	t39 = qJD(3) * t36;
	t38 = qJD(3) * t37;
	t33 = -t35 * t40 + t36 * t41;
	t32 = t34 * t41 + t35 * t39;
	t31 = t34 * t38 + t36 * t42;
	t30 = -t34 * t42 + t36 * t38;
	t1 = [t30, 0, t33, 0, 0, 0; t32, 0, t31, 0, 0, 0; 0, 0, -t39, 0, 0, 0; -t31, 0, -t32, 0, 0, 0; t33, 0, t30, 0, 0, 0; 0, 0, t40, 0, 0, 0; -t41, 0, 0, 0, 0, 0; -t42, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:25:42
	% EndTime: 2019-10-10 00:25:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t47 = sin(qJ(1));
	t52 = qJD(1) * t47;
	t48 = cos(qJ(1));
	t51 = qJD(1) * t48;
	t50 = qJD(3) * t47;
	t49 = qJD(3) * t48;
	t46 = qJ(3) + pkin(9);
	t45 = cos(t46);
	t44 = sin(t46);
	t43 = -t44 * t50 + t45 * t51;
	t42 = t44 * t51 + t45 * t50;
	t41 = t44 * t49 + t45 * t52;
	t40 = -t44 * t52 + t45 * t49;
	t1 = [t40, 0, t43, 0, 0, 0; t42, 0, t41, 0, 0, 0; 0, 0, -qJD(3) * t45, 0, 0, 0; -t41, 0, -t42, 0, 0, 0; t43, 0, t40, 0, 0, 0; 0, 0, qJD(3) * t44, 0, 0, 0; -t51, 0, 0, 0, 0, 0; -t52, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:25:43
	% EndTime: 2019-10-10 00:25:43
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t179 = sin(qJ(1));
	t184 = qJD(1) * t179;
	t180 = cos(qJ(1));
	t183 = qJD(1) * t180;
	t182 = qJD(3) * t179;
	t181 = qJD(3) * t180;
	t178 = qJ(3) + pkin(9);
	t177 = cos(t178);
	t176 = sin(t178);
	t175 = t176 * t182 - t177 * t183;
	t174 = t176 * t183 + t177 * t182;
	t173 = t176 * t181 + t177 * t184;
	t172 = t176 * t184 - t177 * t181;
	t1 = [-t183, 0, 0, 0, 0, 0; -t184, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t172, 0, t175, 0, 0, 0; -t174, 0, -t173, 0, 0, 0; 0, 0, qJD(3) * t177, 0, 0, 0; t173, 0, t174, 0, 0, 0; t175, 0, t172, 0, 0, 0; 0, 0, -qJD(3) * t176, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:25:44
	% EndTime: 2019-10-10 00:25:44
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (102->29), mult. (173->57), div. (0->0), fcn. (173->6), ass. (0->31)
	t260 = cos(qJ(6));
	t261 = cos(qJ(1));
	t280 = t260 * t261;
	t259 = sin(qJ(1));
	t279 = qJD(1) * t259;
	t278 = qJD(1) * t261;
	t258 = sin(qJ(6));
	t277 = qJD(3) * t258;
	t276 = qJD(3) * t259;
	t275 = qJD(3) * t260;
	t274 = qJD(3) * t261;
	t273 = qJD(6) * t258;
	t272 = qJD(6) * t260;
	t271 = qJD(6) * t261;
	t257 = qJ(3) + pkin(9);
	t255 = sin(t257);
	t270 = t255 * t275;
	t269 = t255 * t276;
	t256 = cos(t257);
	t268 = t256 * t276;
	t267 = t255 * t274;
	t266 = t256 * t274;
	t265 = qJD(6) * t256 + qJD(1);
	t264 = qJD(1) * t256 + qJD(6);
	t263 = t265 * t258;
	t262 = t264 * t259 + t267;
	t254 = t262 * t258 - t265 * t280;
	t253 = t262 * t260 + t261 * t263;
	t252 = t265 * t260 * t259 + (t264 * t261 - t269) * t258;
	t251 = -t264 * t280 + (t263 + t270) * t259;
	t1 = [t254, 0, t258 * t268 + (t258 * t278 + t259 * t272) * t255, 0, 0, t251; -t252, 0, -t258 * t266 + (t258 * t279 - t260 * t271) * t255, 0, 0, -t253; 0, 0, -t255 * t277 + t256 * t272, 0, 0, -t255 * t273 + t256 * t275; t253, 0, t260 * t268 + (-t259 * t273 + t260 * t278) * t255, 0, 0, t252; t251, 0, -t260 * t266 + (t258 * t271 + t260 * t279) * t255, 0, 0, t254; 0, 0, -t256 * t273 - t270, 0, 0, -t255 * t272 - t256 * t277; -t255 * t279 + t266, 0, t256 * t278 - t269, 0, 0, 0; t255 * t278 + t268, 0, t256 * t279 + t267, 0, 0, 0; 0, 0, -qJD(3) * t256, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end