% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPRP6
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 16:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPPRP6_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRP6_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_jacobiRD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:07:55
	% EndTime: 2019-12-29 16:07:55
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:07:50
	% EndTime: 2019-12-29 16:07:50
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
	% StartTime: 2019-12-29 16:07:49
	% EndTime: 2019-12-29 16:07:50
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t13 = qJD(1) * sin(qJ(1));
	t11 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; t11, 0, 0, 0, 0; t13, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t13, 0, 0, 0, 0; t11, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:07:50
	% EndTime: 2019-12-29 16:07:50
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (4->4), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t24 = qJD(1) * sin(qJ(1));
	t23 = qJD(1) * cos(qJ(1));
	t20 = cos(pkin(7));
	t19 = sin(pkin(7));
	t1 = [-t19 * t24, 0, 0, 0, 0; t19 * t23, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t20 * t24, 0, 0, 0, 0; t20 * t23, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t23, 0, 0, 0, 0; -t24, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:07:50
	% EndTime: 2019-12-29 16:07:50
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t44 = sin(qJ(1));
	t49 = qJD(1) * t44;
	t45 = cos(qJ(1));
	t48 = qJD(1) * t45;
	t47 = qJD(4) * t44;
	t46 = qJD(4) * t45;
	t43 = pkin(7) + qJ(4);
	t42 = cos(t43);
	t41 = sin(t43);
	t40 = -t41 * t47 + t42 * t48;
	t39 = t41 * t48 + t42 * t47;
	t38 = t41 * t46 + t42 * t49;
	t37 = -t41 * t49 + t42 * t46;
	t1 = [t37, 0, 0, t40, 0; t39, 0, 0, t38, 0; 0, 0, 0, -qJD(4) * t42, 0; -t38, 0, 0, -t39, 0; t40, 0, 0, t37, 0; 0, 0, 0, qJD(4) * t41, 0; -t48, 0, 0, 0, 0; -t49, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:07:51
	% EndTime: 2019-12-29 16:07:51
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (30->11), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t169 = sin(qJ(1));
	t174 = qJD(1) * t169;
	t170 = cos(qJ(1));
	t173 = qJD(1) * t170;
	t172 = qJD(4) * t169;
	t171 = qJD(4) * t170;
	t168 = pkin(7) + qJ(4);
	t167 = cos(t168);
	t166 = sin(t168);
	t165 = -t166 * t172 + t167 * t173;
	t164 = t166 * t173 + t167 * t172;
	t163 = t166 * t171 + t167 * t174;
	t162 = t166 * t174 - t167 * t171;
	t1 = [-t162, 0, 0, t165, 0; t164, 0, 0, t163, 0; 0, 0, 0, -qJD(4) * t167, 0; -t173, 0, 0, 0, 0; -t174, 0, 0, 0, 0; 0, 0, 0, 0, 0; t163, 0, 0, t164, 0; -t165, 0, 0, t162, 0; 0, 0, 0, -qJD(4) * t166, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end