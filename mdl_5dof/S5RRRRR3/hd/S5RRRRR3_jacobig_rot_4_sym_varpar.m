% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
%
% Output:
% Jg_rot [3x5]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S5RRRRR3_jacobig_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_jacobig_rot_4_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_jacobig_rot_4_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:19:38
% EndTime: 2019-07-18 17:19:38
% DurationCPUTime: 0.02s
% Computational Cost: add. (6->4), mult. (2->2), div. (0->0), fcn. (9->4), ass. (0->5)
t87 = cos(qJ(1));
t86 = sin(qJ(1));
t85 = qJ(2) + qJ(3);
t84 = sin(t85);
t1 = [0, t86, t86, t87 * t84, 0; 0, -t87, -t87, t86 * t84, 0; 1, 0, 0, -cos(t85), 0;];
Jg_rot  = t1;
