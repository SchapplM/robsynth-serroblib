% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRPR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:37:01
	% EndTime: 2020-01-03 11:37:01
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:37:01
	% EndTime: 2020-01-03 11:37:01
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t10 = qJD(1) * sin(qJ(1));
	t9 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0; -t10, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t10, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:37:01
	% EndTime: 2020-01-03 11:37:01
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = qJ(1) + pkin(8);
	t14 = qJD(1) * sin(t12);
	t13 = qJD(1) * cos(t12);
	t1 = [0, 0, 0, 0, 0; -t14, 0, 0, 0, 0; t13, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t13, 0, 0, 0, 0; -t14, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:37:01
	% EndTime: 2020-01-03 11:37:02
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (30->8), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t24 = qJ(1) + pkin(8) + qJ(3);
	t25 = qJD(1) + qJD(3);
	t26 = t25 * sin(t24);
	t21 = t25 * cos(t24);
	t1 = [0, 0, 0, 0, 0; -t26, 0, -t26, 0, 0; t21, 0, t21, 0, 0; 0, 0, 0, 0, 0; -t21, 0, -t21, 0, 0; -t26, 0, -t26, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:37:02
	% EndTime: 2020-01-03 11:37:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (40->6), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->13)
	t85 = qJD(1) + qJD(3);
	t91 = t85 * sin(pkin(9));
	t90 = t85 * cos(pkin(9));
	t84 = qJ(1) + pkin(8) + qJ(3);
	t82 = sin(t84);
	t89 = t82 * t90;
	t83 = cos(t84);
	t88 = t83 * t91;
	t81 = t85 * t83;
	t80 = t85 * t82;
	t79 = t83 * t90;
	t78 = t82 * t91;
	t1 = [0, 0, 0, 0, 0; -t89, 0, -t89, 0, 0; t79, 0, t79, 0, 0; 0, 0, 0, 0, 0; t78, 0, t78, 0, 0; -t88, 0, -t88, 0, 0; 0, 0, 0, 0, 0; t81, 0, t81, 0, 0; t80, 0, t80, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:37:02
	% EndTime: 2020-01-03 11:37:02
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (171->15), mult. (132->24), div. (0->0), fcn. (132->6), ass. (0->23)
	t189 = qJD(1) + qJD(3);
	t190 = sin(pkin(9));
	t202 = t189 * t190;
	t191 = cos(pkin(9));
	t192 = sin(qJ(5));
	t201 = t191 * t192;
	t193 = cos(qJ(5));
	t200 = t191 * t193;
	t199 = qJD(5) * t190;
	t188 = qJ(1) + pkin(8) + qJ(3);
	t186 = sin(t188);
	t198 = t186 * t202;
	t187 = cos(t188);
	t197 = t186 * t192 + t187 * t200;
	t196 = t186 * t193 - t187 * t201;
	t195 = -t186 * t200 + t187 * t192;
	t194 = t186 * t201 + t187 * t193;
	t185 = t187 * t202;
	t184 = -t194 * qJD(5) + t197 * t189;
	t183 = t195 * qJD(5) + t196 * t189;
	t182 = t196 * qJD(5) + t195 * t189;
	t181 = -t197 * qJD(5) + t194 * t189;
	t1 = [0, 0, 0, 0, -t193 * t199; t182, 0, t182, 0, t183; t184, 0, t184, 0, -t181; 0, 0, 0, 0, t192 * t199; t181, 0, t181, 0, -t184; t183, 0, t183, 0, t182; 0, 0, 0, 0, 0; -t198, 0, -t198, 0, 0; t185, 0, t185, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end