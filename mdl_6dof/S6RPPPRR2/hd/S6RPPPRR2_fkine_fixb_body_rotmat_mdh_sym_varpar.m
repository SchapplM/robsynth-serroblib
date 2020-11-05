% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPPRR2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:23
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPPPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:23:12
	% EndTime: 2020-11-04 21:23:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:23:12
	% EndTime: 2020-11-04 21:23:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t53 = cos(qJ(1));
	t52 = sin(qJ(1));
	t1 = [t53, -t52, 0, 0; t52, t53, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:23:12
	% EndTime: 2020-11-04 21:23:12
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t56 = qJ(1) + pkin(9);
	t55 = cos(t56);
	t54 = sin(t56);
	t1 = [t55, -t54, 0, cos(qJ(1)) * pkin(1) + 0; t54, t55, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:23:12
	% EndTime: 2020-11-04 21:23:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (6->6), div. (0->0), fcn. (10->4), ass. (0->4)
	t59 = qJ(1) + pkin(9);
	t58 = cos(t59);
	t57 = sin(t59);
	t1 = [0, -t58, t57, t58 * pkin(2) + t57 * qJ(3) + cos(qJ(1)) * pkin(1) + 0; 0, -t57, -t58, t57 * pkin(2) - t58 * qJ(3) + sin(qJ(1)) * pkin(1) + 0; 1, 0, 0, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:23:12
	% EndTime: 2020-11-04 21:23:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (24->14), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->7)
	t65 = pkin(2) + qJ(4);
	t64 = cos(pkin(10));
	t63 = sin(pkin(10));
	t62 = qJ(1) + pkin(9);
	t61 = cos(t62);
	t60 = sin(t62);
	t1 = [t60 * t63, t60 * t64, t61, t65 * t61 + t60 * qJ(3) + cos(qJ(1)) * pkin(1) + 0; -t61 * t63, -t61 * t64, t60, t65 * t60 - t61 * qJ(3) + sin(qJ(1)) * pkin(1) + 0; t64, -t63, 0, pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:23:12
	% EndTime: 2020-11-04 21:23:12
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (39->17), mult. (17->12), div. (0->0), fcn. (25->8), ass. (0->9)
	t75 = pkin(2) + pkin(7) + qJ(4);
	t74 = sin(pkin(10)) * pkin(4) + qJ(3);
	t71 = qJ(1) + pkin(9);
	t70 = pkin(10) + qJ(5);
	t69 = cos(t71);
	t68 = cos(t70);
	t67 = sin(t71);
	t66 = sin(t70);
	t1 = [t67 * t66, t67 * t68, t69, cos(qJ(1)) * pkin(1) + 0 + t75 * t69 + t74 * t67; -t69 * t66, -t69 * t68, t67, sin(qJ(1)) * pkin(1) + 0 - t74 * t69 + t75 * t67; t68, -t66, 0, cos(pkin(10)) * pkin(4) + pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:23:12
	% EndTime: 2020-11-04 21:23:12
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (65->24), mult. (39->24), div. (0->0), fcn. (52->10), ass. (0->15)
	t91 = pkin(2) + pkin(7) + qJ(4);
	t81 = qJ(1) + pkin(9);
	t77 = sin(t81);
	t84 = sin(qJ(6));
	t90 = t77 * t84;
	t85 = cos(qJ(6));
	t89 = t77 * t85;
	t79 = cos(t81);
	t88 = t79 * t84;
	t87 = t79 * t85;
	t80 = pkin(10) + qJ(5);
	t76 = sin(t80);
	t78 = cos(t80);
	t86 = sin(pkin(10)) * pkin(4) + pkin(5) * t76 - pkin(8) * t78 + qJ(3);
	t1 = [t76 * t89 + t88, -t76 * t90 + t87, -t77 * t78, cos(qJ(1)) * pkin(1) + 0 + t91 * t79 + t86 * t77; -t76 * t87 + t90, t76 * t88 + t89, t79 * t78, sin(qJ(1)) * pkin(1) + 0 + t91 * t77 - t86 * t79; t78 * t85, -t78 * t84, t76, t78 * pkin(5) + t76 * pkin(8) + cos(pkin(10)) * pkin(4) + pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end