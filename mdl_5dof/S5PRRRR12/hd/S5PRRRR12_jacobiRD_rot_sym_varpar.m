% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRRR12
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
%   Siehe auch: S5PRRRR12_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PRRRR12_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR12_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
JRD_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t61 = cos(pkin(5));
	t62 = sin(qJ(2));
	t66 = t61 * t62;
	t63 = cos(qJ(2));
	t65 = t61 * t63;
	t64 = qJD(2) * sin(pkin(5));
	t60 = cos(pkin(11));
	t58 = sin(pkin(11));
	t1 = [0, (t58 * t66 - t60 * t63) * qJD(2), 0, 0, 0; 0, (-t58 * t63 - t60 * t66) * qJD(2), 0, 0, 0; 0, -t62 * t64, 0, 0, 0; 0, (t58 * t65 + t60 * t62) * qJD(2), 0, 0, 0; 0, (t58 * t62 - t60 * t65) * qJD(2), 0, 0, 0; 0, -t63 * t64, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (108->12), mult. (72->19), div. (0->0), fcn. (48->8), ass. (0->24)
	t152 = qJD(2) + qJD(3);
	t160 = -t152 / 0.2e1;
	t153 = qJ(2) + qJ(3);
	t149 = pkin(5) - t153;
	t159 = t152 * sin(t149);
	t148 = pkin(5) + t153;
	t158 = t152 * cos(t148);
	t154 = sin(pkin(11));
	t157 = t152 * t154;
	t155 = cos(pkin(11));
	t156 = t152 * t155;
	t151 = cos(t153);
	t150 = sin(t153);
	t145 = cos(t149) * t160;
	t144 = sin(t148) * t160;
	t143 = t145 - t158 / 0.2e1;
	t142 = t159 / 0.2e1 + t144;
	t141 = t158 / 0.2e1 + t145;
	t140 = t144 - t159 / 0.2e1;
	t139 = -t154 * t142 - t151 * t156;
	t138 = -t154 * t143 + t150 * t156;
	t137 = t155 * t142 - t151 * t157;
	t136 = t155 * t143 + t150 * t157;
	t1 = [0, t139, t139, 0, 0; 0, t137, t137, 0, 0; 0, t141, t141, 0, 0; 0, t138, t138, 0, 0; 0, t136, t136, 0, 0; 0, t140, t140, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (258->12), mult. (108->19), div. (0->0), fcn. (72->8), ass. (0->24)
	t178 = qJD(2) + qJD(3) + qJD(4);
	t186 = -t178 / 0.2e1;
	t179 = qJ(2) + qJ(3) + qJ(4);
	t175 = pkin(5) - t179;
	t185 = t178 * sin(t175);
	t174 = pkin(5) + t179;
	t184 = t178 * cos(t174);
	t180 = sin(pkin(11));
	t183 = t178 * t180;
	t181 = cos(pkin(11));
	t182 = t178 * t181;
	t177 = cos(t179);
	t176 = sin(t179);
	t171 = cos(t175) * t186;
	t170 = sin(t174) * t186;
	t169 = t171 - t184 / 0.2e1;
	t168 = t185 / 0.2e1 + t170;
	t167 = t184 / 0.2e1 + t171;
	t166 = t170 - t185 / 0.2e1;
	t165 = -t180 * t168 - t177 * t182;
	t164 = -t180 * t169 + t176 * t182;
	t163 = t181 * t168 - t177 * t183;
	t162 = t181 * t169 + t176 * t183;
	t1 = [0, t165, t165, t165, 0; 0, t163, t163, t163, 0; 0, t167, t167, t167, 0; 0, t164, t164, t164, 0; 0, t162, t162, t162, 0; 0, t166, t166, t166, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:14
	% EndTime: 2024-09-28 18:09:14
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (496->42), mult. (601->94), div. (0->0), fcn. (633->10), ass. (0->55)
	t552 = qJ(2) + qJ(3) + qJ(4);
	t550 = cos(t552);
	t558 = cos(pkin(5));
	t587 = t550 * t558;
	t551 = qJD(2) + qJD(3) + qJD(4);
	t554 = sin(pkin(6));
	t586 = t551 * t554;
	t557 = cos(pkin(6));
	t559 = sin(qJ(5));
	t585 = t557 * t559;
	t560 = cos(qJ(5));
	t584 = t557 * t560;
	t583 = t558 * t559;
	t582 = t558 * t560;
	t581 = qJD(5) * t554;
	t580 = t557 * t583;
	t579 = t557 * t582;
	t578 = t559 * t581;
	t577 = t560 * t581;
	t553 = sin(pkin(11));
	t556 = cos(pkin(11));
	t540 = t553 * t580 - t556 * t560;
	t547 = t553 * t583 - t556 * t584;
	t576 = t540 * qJD(5) + t547 * t551;
	t546 = t553 * t582 + t556 * t585;
	t561 = t553 * t579 + t556 * t559;
	t575 = t561 * qJD(5) + t546 * t551;
	t542 = t553 * t560 + t556 * t580;
	t562 = t553 * t584 + t556 * t583;
	t574 = -t542 * qJD(5) - t551 * t562;
	t543 = -t553 * t559 + t556 * t579;
	t544 = t553 * t585 - t556 * t582;
	t573 = -t543 * qJD(5) + t544 * t551;
	t572 = t544 * qJD(5) - t543 * t551;
	t571 = t562 * qJD(5) + t542 * t551;
	t570 = t546 * qJD(5) + t551 * t561;
	t569 = t547 * qJD(5) + t540 * t551;
	t555 = sin(pkin(5));
	t568 = t555 * t578;
	t567 = t555 * t577;
	t549 = sin(t552);
	t566 = t549 * t559 - t550 * t584;
	t565 = -t549 * t560 - t550 * t585;
	t564 = -t549 * t584 - t550 * t559;
	t563 = t549 * t585 - t550 * t560;
	t548 = t555 * t550 * t586;
	t531 = (-t549 * t553 + t556 * t587) * t586;
	t530 = (-t549 * t556 - t553 * t587) * t586;
	t529 = (t564 * qJD(5) + t565 * t551) * t555;
	t528 = (t563 * qJD(5) + t566 * t551) * t555;
	t527 = t573 * t549 - t571 * t550;
	t526 = -t574 * t549 + t572 * t550;
	t525 = t575 * t549 + t569 * t550;
	t524 = -t576 * t549 + t570 * t550;
	t1 = [0, t525, t525, t525, t570 * t549 + t576 * t550 - t553 * t568; 0, t527, t527, t527, t572 * t549 + t574 * t550 + t556 * t568; 0, t529, t529, t529, -t558 * t578 + (t565 * qJD(5) + t564 * t551) * t555; 0, t524, t524, t524, -t569 * t549 + t575 * t550 - t553 * t567; 0, t526, t526, t526, t571 * t549 + t573 * t550 + t556 * t567; 0, t528, t528, t528, -t558 * t577 + (t566 * qJD(5) + t563 * t551) * t555; 0, t530, t530, t530, 0; 0, t531, t531, t531, 0; 0, t548, t548, t548, 0;];
	JRD_rot = t1;
end