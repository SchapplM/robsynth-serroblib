% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR11 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:21
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRPR11_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR11_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:21:00
	% EndTime: 2020-11-04 20:21:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:21:00
	% EndTime: 2020-11-04 20:21:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t44 = cos(qJ(1));
	t43 = sin(qJ(1));
	t1 = [t44, -t43, 0, 0; t43, t44, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:21:00
	% EndTime: 2020-11-04 20:21:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t48 = cos(qJ(1));
	t47 = sin(qJ(1));
	t46 = cos(pkin(8));
	t45 = sin(pkin(8));
	t1 = [t48 * t46, -t48 * t45, t47, t48 * pkin(1) + t47 * qJ(2) + 0; t47 * t46, -t47 * t45, -t48, t47 * pkin(1) - t48 * qJ(2) + 0; t45, t46, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:21:00
	% EndTime: 2020-11-04 20:21:00
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t55 = cos(qJ(1));
	t54 = sin(qJ(1));
	t53 = pkin(6) + qJ(2);
	t52 = pkin(8) + qJ(3);
	t51 = cos(t52);
	t50 = sin(t52);
	t49 = cos(pkin(8)) * pkin(2) + pkin(1);
	t1 = [t55 * t51, -t55 * t50, t54, t49 * t55 + t53 * t54 + 0; t54 * t51, -t54 * t50, -t55, t49 * t54 - t53 * t55 + 0; t50, t51, 0, sin(pkin(8)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:21:00
	% EndTime: 2020-11-04 20:21:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (30->15), mult. (21->14), div. (0->0), fcn. (29->6), ass. (0->8)
	t59 = pkin(8) + qJ(3);
	t57 = sin(t59);
	t58 = cos(t59);
	t63 = pkin(3) * t58 + qJ(4) * t57 + cos(pkin(8)) * pkin(2) + pkin(1);
	t62 = cos(qJ(1));
	t61 = sin(qJ(1));
	t60 = pkin(6) + qJ(2);
	t1 = [t62 * t58, t61, t62 * t57, t60 * t61 + t63 * t62 + 0; t61 * t58, -t62, t61 * t57, -t62 * t60 + t63 * t61 + 0; t57, 0, -t58, t57 * pkin(3) - t58 * qJ(4) + sin(pkin(8)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:21:00
	% EndTime: 2020-11-04 20:21:00
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (48->21), mult. (37->22), div. (0->0), fcn. (51->10), ass. (0->15)
	t79 = cos(qJ(5));
	t78 = sin(qJ(5));
	t77 = pkin(3) + pkin(4);
	t76 = cos(qJ(1));
	t75 = sin(qJ(1));
	t74 = cos(pkin(8));
	t73 = sin(pkin(8));
	t72 = pkin(8) + qJ(3);
	t71 = qJ(2) + pkin(6) - pkin(7);
	t70 = cos(t72);
	t69 = sin(t72);
	t66 = t69 * t79 - t70 * t78;
	t65 = -t69 * t78 - t70 * t79;
	t64 = (qJ(4) * t73 + t77 * t74) * cos(qJ(3)) + (qJ(4) * t74 - t73 * t77) * sin(qJ(3)) + t74 * pkin(2) + pkin(1);
	t1 = [-t76 * t65, t76 * t66, -t75, t64 * t76 + t71 * t75 + 0; -t75 * t65, t75 * t66, t76, t64 * t75 - t71 * t76 + 0; t66, t65, 0, t73 * pkin(2) - t70 * qJ(4) + t77 * t69 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end