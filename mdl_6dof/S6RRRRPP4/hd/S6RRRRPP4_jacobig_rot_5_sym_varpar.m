% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP4
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRPP4_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_jacobig_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:26:59
% EndTime: 2019-02-26 22:26:59
% DurationCPUTime: 0.02s
% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (12->4), ass. (0->7)
t78 = cos(qJ(1));
t77 = cos(qJ(2));
t76 = sin(qJ(1));
t75 = sin(qJ(2));
t74 = t78 * t75;
t73 = t76 * t75;
t1 = [0, t76, t74, t74, 0, 0; 0, -t78, t73, t73, 0, 0; 1, 0, -t77, -t77, 0, 0;];
Jg_rot  = t1;