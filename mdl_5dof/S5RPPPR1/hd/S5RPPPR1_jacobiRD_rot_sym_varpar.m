% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:38
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPPPR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:38:58
	% EndTime: 2019-10-24 10:38:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:38:58
	% EndTime: 2019-10-24 10:38:58
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = qJD(1) * cos(qJ(1));
	t7 = qJD(1) * sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; t7, 0, 0, 0, 0; -t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; t9, 0, 0, 0, 0; t7, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:38:58
	% EndTime: 2019-10-24 10:38:58
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (5->2), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = qJ(1) + pkin(7);
	t13 = qJD(1) * cos(t12);
	t10 = qJD(1) * sin(t12);
	t1 = [0, 0, 0, 0, 0; t10, 0, 0, 0, 0; -t13, 0, 0, 0, 0; 0, 0, 0, 0, 0; t13, 0, 0, 0, 0; t10, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:38:58
	% EndTime: 2019-10-24 10:38:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->5), mult. (10->8), div. (0->0), fcn. (10->4), ass. (0->6)
	t49 = qJD(1) * sin(pkin(8));
	t48 = qJD(1) * cos(pkin(8));
	t45 = qJ(1) + pkin(7);
	t44 = cos(t45);
	t43 = sin(t45);
	t1 = [0, 0, 0, 0, 0; t43 * t48, 0, 0, 0, 0; -t44 * t48, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t43 * t49, 0, 0, 0, 0; t44 * t49, 0, 0, 0, 0; 0, 0, 0, 0, 0; -qJD(1) * t44, 0, 0, 0, 0; -qJD(1) * t43, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:38:58
	% EndTime: 2019-10-24 10:38:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (15->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->10)
	t93 = sin(pkin(9));
	t96 = cos(pkin(8));
	t99 = t93 * t96;
	t95 = cos(pkin(9));
	t98 = t95 * t96;
	t97 = qJD(1) * sin(pkin(8));
	t92 = qJ(1) + pkin(7);
	t91 = cos(t92);
	t90 = sin(t92);
	t1 = [0, 0, 0, 0, 0; (t90 * t98 - t91 * t93) * qJD(1), 0, 0, 0, 0; (-t90 * t93 - t91 * t98) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0; (-t90 * t99 - t91 * t95) * qJD(1), 0, 0, 0, 0; (-t90 * t95 + t91 * t99) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0; t90 * t97, 0, 0, 0, 0; -t91 * t97, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:38:58
	% EndTime: 2019-10-24 10:38:58
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (94->14), mult. (88->24), div. (0->0), fcn. (88->6), ass. (0->21)
	t156 = qJ(1) + pkin(7);
	t152 = sin(t156);
	t158 = cos(pkin(8));
	t166 = t152 * t158;
	t154 = cos(t156);
	t165 = t154 * t158;
	t157 = sin(pkin(8));
	t164 = qJD(1) * t157;
	t163 = qJD(5) * t157;
	t155 = pkin(9) + qJ(5);
	t151 = sin(t155);
	t153 = cos(t155);
	t162 = t151 * t152 + t153 * t165;
	t161 = -t151 * t154 + t153 * t166;
	t160 = t151 * t165 - t152 * t153;
	t159 = t151 * t166 + t153 * t154;
	t150 = t162 * qJD(1) - t159 * qJD(5);
	t149 = t160 * qJD(1) + t161 * qJD(5);
	t148 = t161 * qJD(1) + t160 * qJD(5);
	t147 = t159 * qJD(1) - t162 * qJD(5);
	t1 = [0, 0, 0, 0, -t153 * t163; t148, 0, 0, 0, t149; -t150, 0, 0, 0, t147; 0, 0, 0, 0, t151 * t163; -t147, 0, 0, 0, t150; t149, 0, 0, 0, t148; 0, 0, 0, 0, 0; t152 * t164, 0, 0, 0, 0; -t154 * t164, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end