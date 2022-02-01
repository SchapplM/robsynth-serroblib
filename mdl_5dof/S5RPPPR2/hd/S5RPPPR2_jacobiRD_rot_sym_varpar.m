% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPPR2
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
%   Siehe auch: S5RPPPR2_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPPPR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
JRD_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t22 = qJD(1) * sin(qJ(1));
	t21 = qJD(1) * cos(qJ(1));
	t18 = cos(pkin(7));
	t17 = sin(pkin(7));
	t1 = [-t18 * t21, 0, 0, 0, 0; -t18 * t22, 0, 0, 0, 0; 0, 0, 0, 0, 0; t17 * t21, 0, 0, 0, 0; t17 * t22, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t22, 0, 0, 0, 0; t21, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:32
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t128 = cos(pkin(7));
	t129 = sin(qJ(1));
	t133 = t128 * t129;
	t130 = cos(qJ(1));
	t132 = t128 * t130;
	t131 = qJD(1) * sin(pkin(7));
	t127 = cos(pkin(8));
	t125 = sin(pkin(8));
	t1 = [(-t125 * t129 - t127 * t132) * qJD(1), 0, 0, 0, 0; (t125 * t130 - t127 * t133) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0; (t125 * t132 - t127 * t129) * qJD(1), 0, 0, 0, 0; (t125 * t133 + t127 * t130) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0; -t130 * t131, 0, 0, 0, 0; -t129 * t131, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:32
	% EndTime: 2022-01-23 09:00:32
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (10->8), mult. (46->25), div. (0->0), fcn. (46->8), ass. (0->15)
	t188 = cos(pkin(7));
	t189 = sin(qJ(1));
	t195 = t188 * t189;
	t190 = cos(qJ(1));
	t194 = t188 * t190;
	t193 = qJD(1) * sin(pkin(7));
	t192 = t189 * t193;
	t191 = t190 * t193;
	t187 = cos(pkin(8));
	t186 = cos(pkin(9));
	t184 = sin(pkin(8));
	t183 = sin(pkin(9));
	t182 = (-t184 * t189 - t187 * t194) * qJD(1);
	t181 = (t184 * t190 - t187 * t195) * qJD(1);
	t1 = [t182 * t186 - t183 * t191, 0, 0, 0, 0; t181 * t186 - t183 * t192, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t182 * t183 - t186 * t191, 0, 0, 0, 0; -t181 * t183 - t186 * t192, 0, 0, 0, 0; 0, 0, 0, 0, 0; (-t184 * t194 + t187 * t189) * qJD(1), 0, 0, 0, 0; (-t184 * t195 - t187 * t190) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:32
	% EndTime: 2022-01-23 09:00:32
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (80->32), mult. (258->74), div. (0->0), fcn. (292->10), ass. (0->33)
	t332 = sin(pkin(9));
	t334 = sin(pkin(7));
	t337 = cos(pkin(7));
	t335 = cos(pkin(9));
	t336 = cos(pkin(8));
	t354 = t335 * t336;
	t328 = t334 * t332 + t337 * t354;
	t338 = sin(qJ(5));
	t333 = sin(pkin(8));
	t340 = cos(qJ(5));
	t356 = t333 * t340;
	t363 = -t328 * t338 + t337 * t356;
	t320 = t363 * qJD(5);
	t357 = t333 * t338;
	t329 = t335 * t357 + t340 * t336;
	t325 = t329 * qJD(5);
	t330 = t335 * t356 - t338 * t336;
	t339 = sin(qJ(1));
	t341 = cos(qJ(1));
	t344 = t328 * t340 + t337 * t357;
	t365 = -t341 * t325 - t320 * t339 + (-t339 * t330 - t344 * t341) * qJD(1);
	t347 = t341 * t336;
	t351 = t337 * t339;
	t352 = t335 * t341;
	t362 = -(t333 * t351 + t347) * t338 + (-t328 * t339 + t333 * t352) * t340;
	t321 = t344 * qJD(5);
	t326 = t330 * qJD(5);
	t361 = -t321 * t341 - t339 * t326 + (-t341 * t329 - t339 * t363) * qJD(1);
	t355 = t333 * t341;
	t353 = t335 * t339;
	t342 = -(t328 * t341 + t333 * t353) * t338 + (-t336 * t339 + t337 * t355) * t340;
	t327 = -t337 * t332 + t334 * t354;
	t1 = [t365, 0, 0, 0, t361; t362 * qJD(1) + t342 * qJD(5), 0, 0, 0, t342 * qJD(1) + t362 * qJD(5); 0, 0, 0, 0, (-t327 * t340 - t334 * t357) * qJD(5); t321 * t339 - t341 * t326 + (t329 * t339 - t341 * t363) * qJD(1), 0, 0, 0, -t320 * t341 + t339 * t325 + (-t330 * t341 + t339 * t344) * qJD(1); t361, 0, 0, 0, t365; 0, 0, 0, 0, (t327 * t338 - t334 * t356) * qJD(5); ((-t333 * t339 - t337 * t347) * t332 + t334 * t352) * qJD(1), 0, 0, 0, 0; ((-t336 * t351 + t355) * t332 + t334 * t353) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
end