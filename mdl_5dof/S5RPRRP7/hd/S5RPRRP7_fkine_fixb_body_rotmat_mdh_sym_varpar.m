% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:24
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRRP7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:24:15
	% EndTime: 2020-11-04 20:24:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:24:15
	% EndTime: 2020-11-04 20:24:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t51 = cos(qJ(1));
	t50 = sin(qJ(1));
	t1 = [t51, -t50, 0, 0; t50, t51, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:24:15
	% EndTime: 2020-11-04 20:24:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t54 = qJ(1) + pkin(8);
	t53 = cos(t54);
	t52 = sin(t54);
	t1 = [t53, -t52, 0, cos(qJ(1)) * pkin(1) + 0; t52, t53, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:24:15
	% EndTime: 2020-11-04 20:24:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t59 = cos(qJ(3));
	t58 = sin(qJ(3));
	t57 = qJ(1) + pkin(8);
	t56 = cos(t57);
	t55 = sin(t57);
	t1 = [t56 * t59, -t56 * t58, t55, t56 * pkin(2) + t55 * pkin(6) + cos(qJ(1)) * pkin(1) + 0; t55 * t59, -t55 * t58, -t56, t55 * pkin(2) - t56 * pkin(6) + sin(qJ(1)) * pkin(1) + 0; t58, t59, 0, qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:24:15
	% EndTime: 2020-11-04 20:24:15
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (38->19), mult. (32->24), div. (0->0), fcn. (45->8), ass. (0->11)
	t63 = sin(qJ(4));
	t66 = cos(qJ(3));
	t69 = t63 * t66;
	t65 = cos(qJ(4));
	t68 = t65 * t66;
	t64 = sin(qJ(3));
	t67 = pkin(3) * t66 + pkin(7) * t64 + pkin(2);
	t62 = qJ(1) + pkin(8);
	t61 = cos(t62);
	t60 = sin(t62);
	t1 = [t60 * t63 + t61 * t68, t60 * t65 - t61 * t69, t61 * t64, cos(qJ(1)) * pkin(1) + t60 * pkin(6) + 0 + t67 * t61; t60 * t68 - t61 * t63, -t60 * t69 - t61 * t65, t60 * t64, sin(qJ(1)) * pkin(1) - t61 * pkin(6) + 0 + t67 * t60; t64 * t65, -t64 * t63, -t66, t64 * pkin(3) - t66 * pkin(7) + pkin(5) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:24:15
	% EndTime: 2020-11-04 20:24:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (55->24), mult. (52->30), div. (0->0), fcn. (69->8), ass. (0->15)
	t77 = sin(qJ(4));
	t80 = cos(qJ(3));
	t83 = t77 * t80;
	t79 = cos(qJ(4));
	t82 = t79 * t80;
	t78 = sin(qJ(3));
	t81 = pkin(3) * t80 + pkin(7) * t78 + pkin(2);
	t76 = qJ(1) + pkin(8);
	t75 = cos(t76);
	t74 = sin(t76);
	t73 = t74 * t77 + t75 * t82;
	t72 = -t74 * t79 + t75 * t83;
	t71 = t74 * t82 - t75 * t77;
	t70 = t74 * t83 + t75 * t79;
	t1 = [t73, t75 * t78, t72, cos(qJ(1)) * pkin(1) + t73 * pkin(4) + t74 * pkin(6) + t72 * qJ(5) + 0 + t81 * t75; t71, t74 * t78, t70, sin(qJ(1)) * pkin(1) + t71 * pkin(4) - t75 * pkin(6) + t70 * qJ(5) + 0 + t81 * t74; t78 * t79, -t80, t78 * t77, -t80 * pkin(7) + pkin(5) + qJ(2) + 0 + (pkin(4) * t79 + qJ(5) * t77 + pkin(3)) * t78; 0, 0, 0, 1;];
	Tc_mdh = t1;
end