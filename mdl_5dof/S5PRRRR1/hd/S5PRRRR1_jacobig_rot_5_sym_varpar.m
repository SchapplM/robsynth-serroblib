% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
%
% Output:
% Jg_rot [3x5]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S5PRRRR1_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_jacobig_rot_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_jacobig_rot_5_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:29:14
% EndTime: 2019-07-18 13:29:14
% DurationCPUTime: 0.02s
% Computational Cost: add. (5->3), mult. (2->2), div. (0->0), fcn. (9->4), ass. (0->5)
t81 = cos(qJ(2));
t80 = sin(qJ(2));
t79 = qJ(3) + qJ(4);
t78 = sin(t79);
t1 = [0, 0, t80, t80, t81 * t78; 0, -1, 0, 0, cos(t79); 0, 0, -t81, -t81, t80 * t78;];
Jg_rot  = t1;
