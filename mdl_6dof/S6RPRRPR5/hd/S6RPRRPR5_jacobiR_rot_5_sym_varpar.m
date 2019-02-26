% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR5
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPR5_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:03:12
% EndTime: 2019-02-26 21:03:12
% DurationCPUTime: 0.02s
% Computational Cost: add. (36->5), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
t77 = cos(qJ(1));
t76 = sin(qJ(1));
t75 = pkin(10) + qJ(3) + qJ(4);
t74 = cos(t75);
t73 = sin(t75);
t72 = t77 * t74;
t71 = t77 * t73;
t70 = t76 * t74;
t69 = t76 * t73;
t1 = [t77, 0, 0, 0, 0, 0; t76, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t70, 0, t71, t71, 0, 0; -t72, 0, t69, t69, 0, 0; 0, 0, -t74, -t74, 0, 0; -t69, 0, t72, t72, 0, 0; t71, 0, t70, t70, 0, 0; 0, 0, t73, t73, 0, 0;];
JR_rot  = t1;
