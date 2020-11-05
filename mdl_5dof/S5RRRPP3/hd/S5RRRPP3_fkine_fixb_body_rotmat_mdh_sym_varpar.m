% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRPP3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:40
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRPP3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:40:11
	% EndTime: 2020-11-04 20:40:11
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:40:11
	% EndTime: 2020-11-04 20:40:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t41 = cos(qJ(1));
	t40 = sin(qJ(1));
	t1 = [t41, -t40, 0, 0; t40, t41, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:40:11
	% EndTime: 2020-11-04 20:40:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t44 = qJ(1) + qJ(2);
	t43 = cos(t44);
	t42 = sin(t44);
	t1 = [t43, -t42, 0, cos(qJ(1)) * pkin(1) + 0; t42, t43, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:40:11
	% EndTime: 2020-11-04 20:40:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t49 = cos(qJ(3));
	t48 = sin(qJ(3));
	t47 = qJ(1) + qJ(2);
	t46 = cos(t47);
	t45 = sin(t47);
	t1 = [t46 * t49, -t46 * t48, t45, t46 * pkin(2) + t45 * pkin(7) + cos(qJ(1)) * pkin(1) + 0; t45 * t49, -t45 * t48, -t46, t45 * pkin(2) - t46 * pkin(7) + sin(qJ(1)) * pkin(1) + 0; t48, t49, 0, pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:40:11
	% EndTime: 2020-11-04 20:40:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->18), mult. (20->14), div. (0->0), fcn. (28->6), ass. (0->7)
	t53 = sin(qJ(3));
	t54 = cos(qJ(3));
	t55 = pkin(3) * t54 + qJ(4) * t53 + pkin(2);
	t52 = qJ(1) + qJ(2);
	t51 = cos(t52);
	t50 = sin(t52);
	t1 = [t50, -t51 * t54, t51 * t53, cos(qJ(1)) * pkin(1) + t50 * pkin(7) + 0 + t55 * t51; -t51, -t50 * t54, t50 * t53, sin(qJ(1)) * pkin(1) - t51 * pkin(7) + 0 + t55 * t50; 0, -t53, -t54, t53 * pkin(3) - t54 * qJ(4) + pkin(5) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:40:11
	% EndTime: 2020-11-04 20:40:11
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (39->17), mult. (27->14), div. (0->0), fcn. (35->6), ass. (0->9)
	t63 = pkin(4) + pkin(7);
	t62 = pkin(3) + qJ(5);
	t59 = sin(qJ(3));
	t60 = cos(qJ(3));
	t61 = qJ(4) * t59 + t62 * t60 + pkin(2);
	t58 = qJ(1) + qJ(2);
	t57 = cos(t58);
	t56 = sin(t58);
	t1 = [t56, t57 * t59, t57 * t60, cos(qJ(1)) * pkin(1) + 0 + t63 * t56 + t61 * t57; -t57, t56 * t59, t56 * t60, sin(qJ(1)) * pkin(1) + 0 - t63 * t57 + t61 * t56; 0, -t60, t59, -t60 * qJ(4) + t62 * t59 + pkin(5) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end