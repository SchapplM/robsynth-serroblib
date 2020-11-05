% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PPRRR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:56
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PPRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:56:09
	% EndTime: 2020-11-04 19:56:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:56:09
	% EndTime: 2020-11-04 19:56:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t42 = cos(pkin(8));
	t41 = sin(pkin(8));
	t1 = [t42, -t41, 0, 0; t41, t42, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:56:09
	% EndTime: 2020-11-04 19:56:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (6->6), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t44 = cos(pkin(8));
	t43 = sin(pkin(8));
	t1 = [t44, 0, t43, t44 * pkin(1) + t43 * qJ(2) + 0; t43, 0, -t44, t43 * pkin(1) - t44 * qJ(2) + 0; 0, 1, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:56:09
	% EndTime: 2020-11-04 19:56:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->10), mult. (12->8), div. (0->0), fcn. (20->4), ass. (0->8)
	t51 = pkin(1) + pkin(2);
	t50 = cos(qJ(3));
	t49 = sin(qJ(3));
	t48 = cos(pkin(8));
	t47 = sin(pkin(8));
	t46 = t47 * t50 - t48 * t49;
	t45 = -t47 * t49 - t48 * t50;
	t1 = [-t45, t46, 0, t47 * qJ(2) + t51 * t48 + 0; t46, t45, 0, -t48 * qJ(2) + t51 * t47 + 0; 0, 0, -1, -pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:56:09
	% EndTime: 2020-11-04 19:56:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (25->16), mult. (20->14), div. (0->0), fcn. (28->6), ass. (0->11)
	t61 = pkin(1) + pkin(2);
	t60 = cos(qJ(3));
	t59 = sin(qJ(3));
	t58 = cos(pkin(8));
	t57 = sin(pkin(8));
	t56 = qJ(3) + qJ(4);
	t55 = cos(t56);
	t54 = sin(t56);
	t53 = -t58 * t54 + t57 * t55;
	t52 = -t57 * t54 - t58 * t55;
	t1 = [-t52, t53, 0, t57 * qJ(2) + t61 * t58 + 0 + (t57 * t59 + t58 * t60) * pkin(3); t53, t52, 0, -t58 * qJ(2) + t61 * t57 + 0 + (t57 * t60 - t58 * t59) * pkin(3); 0, 0, -1, -pkin(6) - pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:56:09
	% EndTime: 2020-11-04 19:56:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (46->23), mult. (26->18), div. (0->0), fcn. (40->10), ass. (0->15)
	t69 = pkin(8) - qJ(3);
	t75 = pkin(1) + pkin(2);
	t74 = cos(qJ(5));
	t73 = sin(qJ(5));
	t72 = cos(pkin(8));
	t71 = sin(pkin(8));
	t70 = qJ(3) + qJ(4);
	t68 = cos(t70);
	t67 = sin(t70);
	t66 = -qJ(4) + t69;
	t65 = cos(t66);
	t64 = sin(t66);
	t63 = t71 * t67 + t72 * t68;
	t62 = t72 * t67 - t71 * t68;
	t1 = [t63 * t74, -t63 * t73, t62, pkin(4) * t65 - pkin(7) * t64 + pkin(3) * cos(t69) + t75 * t72 + t71 * qJ(2) + 0; -t62 * t74, t62 * t73, t63, pkin(7) * t65 + pkin(4) * t64 + pkin(3) * sin(t69) + t75 * t71 - t72 * qJ(2) + 0; -t73, -t74, 0, -pkin(6) - pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end