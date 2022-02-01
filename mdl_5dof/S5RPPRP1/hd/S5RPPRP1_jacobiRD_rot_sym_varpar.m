% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPRP1
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
%   Siehe auch: S5RPPRP1_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPPRP1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRP1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
JRD_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:04
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:04
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:04
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = qJ(1) + pkin(7);
	t37 = qJD(1) * sin(t35);
	t36 = qJD(1) * cos(t35);
	t1 = [-t36, 0, 0, 0, 0; -t37, 0, 0, 0, 0; 0, 0, 0, 0, 0; t37, 0, 0, 0, 0; -t36, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:04
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (9->4), mult. (10->8), div. (0->0), fcn. (10->4), ass. (0->6)
	t26 = qJD(1) * sin(pkin(8));
	t25 = qJD(1) * cos(pkin(8));
	t22 = qJ(1) + pkin(7);
	t21 = cos(t22);
	t20 = sin(t22);
	t1 = [-t21 * t25, 0, 0, 0, 0; -t20 * t25, 0, 0, 0, 0; 0, 0, 0, 0, 0; t21 * t26, 0, 0, 0, 0; t20 * t26, 0, 0, 0, 0; 0, 0, 0, 0, 0; -qJD(1) * t20, 0, 0, 0, 0; qJD(1) * t21, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:05
	% EndTime: 2022-01-23 09:13:05
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (61->14), mult. (88->24), div. (0->0), fcn. (88->6), ass. (0->20)
	t162 = cos(pkin(8));
	t163 = sin(qJ(4));
	t172 = t162 * t163;
	t164 = cos(qJ(4));
	t171 = t162 * t164;
	t161 = sin(pkin(8));
	t170 = qJD(1) * t161;
	t169 = qJD(4) * t161;
	t160 = qJ(1) + pkin(7);
	t158 = sin(t160);
	t159 = cos(t160);
	t168 = -t158 * t163 - t159 * t171;
	t167 = -t158 * t164 + t159 * t172;
	t166 = t158 * t171 - t159 * t163;
	t165 = t158 * t172 + t159 * t164;
	t157 = t168 * qJD(1) + t165 * qJD(4);
	t156 = t167 * qJD(1) + t166 * qJD(4);
	t155 = t166 * qJD(1) + t167 * qJD(4);
	t154 = t165 * qJD(1) + t168 * qJD(4);
	t1 = [t157, 0, 0, t154, 0; -t155, 0, 0, -t156, 0; 0, 0, 0, -t164 * t169, 0; t156, 0, 0, t155, 0; t154, 0, 0, t157, 0; 0, 0, 0, t163 * t169, 0; -t159 * t170, 0, 0, 0, 0; -t158 * t170, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:05
	% EndTime: 2022-01-23 09:13:05
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (61->14), mult. (88->24), div. (0->0), fcn. (88->6), ass. (0->20)
	t173 = cos(pkin(8));
	t174 = sin(qJ(4));
	t183 = t173 * t174;
	t175 = cos(qJ(4));
	t182 = t173 * t175;
	t172 = sin(pkin(8));
	t181 = qJD(1) * t172;
	t180 = qJD(4) * t172;
	t171 = qJ(1) + pkin(7);
	t169 = sin(t171);
	t170 = cos(t171);
	t179 = -t169 * t174 - t170 * t182;
	t178 = -t169 * t175 + t170 * t183;
	t177 = t169 * t182 - t170 * t174;
	t176 = t169 * t183 + t170 * t175;
	t168 = t179 * qJD(1) + t176 * qJD(4);
	t167 = t178 * qJD(1) + t177 * qJD(4);
	t166 = t177 * qJD(1) + t178 * qJD(4);
	t165 = t176 * qJD(1) + t179 * qJD(4);
	t1 = [t168, 0, 0, t165, 0; -t166, 0, 0, -t167, 0; 0, 0, 0, -t175 * t180, 0; t167, 0, 0, t166, 0; t165, 0, 0, t168, 0; 0, 0, 0, t174 * t180, 0; -t170 * t181, 0, 0, 0, 0; -t169 * t181, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
end