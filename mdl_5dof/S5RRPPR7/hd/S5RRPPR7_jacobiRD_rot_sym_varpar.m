% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:23
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRPPR7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:23:04
	% EndTime: 2019-12-29 18:23:04
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:23:09
	% EndTime: 2019-12-29 18:23:09
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:23:05
	% EndTime: 2019-12-29 18:23:05
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t34 = sin(qJ(1));
	t41 = qJD(1) * t34;
	t36 = cos(qJ(1));
	t40 = qJD(1) * t36;
	t33 = sin(qJ(2));
	t39 = qJD(2) * t33;
	t35 = cos(qJ(2));
	t38 = qJD(2) * t35;
	t37 = qJD(2) * t36;
	t32 = t34 * t39 - t35 * t40;
	t31 = t33 * t40 + t34 * t38;
	t30 = t33 * t37 + t35 * t41;
	t29 = t33 * t41 - t35 * t37;
	t1 = [t32, t29, 0, 0, 0; -t30, -t31, 0, 0, 0; 0, -t39, 0, 0, 0; t31, t30, 0, 0, 0; t29, t32, 0, 0, 0; 0, -t38, 0, 0, 0; -t41, 0, 0, 0, 0; t40, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:23:05
	% EndTime: 2019-12-29 18:23:05
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t50 = sin(qJ(1));
	t55 = qJD(1) * t50;
	t51 = cos(qJ(1));
	t54 = qJD(1) * t51;
	t53 = qJD(2) * t50;
	t52 = qJD(2) * t51;
	t49 = qJ(2) + pkin(8);
	t48 = cos(t49);
	t47 = sin(t49);
	t46 = t47 * t53 - t48 * t54;
	t45 = t47 * t54 + t48 * t53;
	t44 = t47 * t52 + t48 * t55;
	t43 = t47 * t55 - t48 * t52;
	t1 = [t46, t43, 0, 0, 0; -t44, -t45, 0, 0, 0; 0, -qJD(2) * t47, 0, 0, 0; t45, t44, 0, 0, 0; t43, t46, 0, 0, 0; 0, -qJD(2) * t48, 0, 0, 0; -t55, 0, 0, 0, 0; t54, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:23:11
	% EndTime: 2019-12-29 18:23:11
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (27->8), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t180 = sin(qJ(1));
	t185 = qJD(1) * t180;
	t181 = cos(qJ(1));
	t184 = qJD(1) * t181;
	t183 = qJD(2) * t180;
	t182 = qJD(2) * t181;
	t179 = qJ(2) + pkin(8);
	t178 = cos(t179);
	t177 = sin(t179);
	t176 = -t177 * t183 + t178 * t184;
	t175 = t177 * t184 + t178 * t183;
	t174 = t177 * t182 + t178 * t185;
	t173 = -t177 * t185 + t178 * t182;
	t1 = [-t185, 0, 0, 0, 0; t184, 0, 0, 0, 0; 0, 0, 0, 0, 0; t176, t173, 0, 0, 0; t174, t175, 0, 0, 0; 0, qJD(2) * t177, 0, 0, 0; -t175, -t174, 0, 0, 0; t173, t176, 0, 0, 0; 0, qJD(2) * t178, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:23:12
	% EndTime: 2019-12-29 18:23:12
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (102->29), mult. (173->56), div. (0->0), fcn. (173->6), ass. (0->32)
	t259 = cos(qJ(5));
	t256 = qJ(2) + pkin(8);
	t254 = sin(t256);
	t264 = qJD(5) * t254 + qJD(1);
	t280 = t259 * t264;
	t257 = sin(qJ(5));
	t279 = t264 * t257;
	t258 = sin(qJ(1));
	t278 = qJD(1) * t258;
	t260 = cos(qJ(1));
	t277 = qJD(1) * t260;
	t276 = qJD(2) * t257;
	t275 = qJD(2) * t258;
	t274 = qJD(2) * t259;
	t273 = qJD(2) * t260;
	t272 = qJD(5) * t257;
	t271 = qJD(5) * t259;
	t270 = qJD(5) * t260;
	t255 = cos(t256);
	t269 = t255 * t274;
	t268 = t254 * t275;
	t267 = t255 * t275;
	t266 = t254 * t273;
	t265 = t255 * t273;
	t263 = -qJD(1) * t254 - qJD(5);
	t262 = t263 * t260;
	t261 = t263 * t258 + t265;
	t253 = t261 * t257 + t260 * t280;
	t252 = t261 * t259 - t260 * t279;
	t251 = -t258 * t280 + (t262 - t267) * t257;
	t250 = t259 * t262 + (-t269 + t279) * t258;
	t1 = [t251, -t257 * t266 + (-t257 * t278 + t259 * t270) * t255, 0, 0, t252; t253, -t257 * t268 + (t257 * t277 + t258 * t271) * t255, 0, 0, -t250; 0, t254 * t271 + t255 * t276, 0, 0, t254 * t274 + t255 * t272; t250, -t259 * t266 + (-t257 * t270 - t259 * t278) * t255, 0, 0, -t253; t252, -t259 * t268 + (-t258 * t272 + t259 * t277) * t255, 0, 0, t251; 0, -t254 * t272 + t269, 0, 0, -t254 * t276 + t255 * t271; -t255 * t277 + t268, t254 * t278 - t265, 0, 0, 0; -t255 * t278 - t266, -t254 * t277 - t267, 0, 0, 0; 0, -qJD(2) * t254, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end